#include "qing_mcost_to_disp.h"

//bilateral filter weight
//float * qx_get_color_weighted_table(float sigma_range, int len) {
//    float * table_color, * color_table_x;
//    table_color = new float[len];
//    color_table_x = &table_color[0];
//    for(int y = 0; y<len; ++y) {
//        *color_table_x++ = exp(-double(y*y)/2*sigma_range*sigma_range);
//    }
//    return table_color;
//}

//can be seperated into two-dimensional
//float * qing_get_range_weighted_table(float sigma_range, int len) {
//    float * table_color, * color_table_x;
//    table_color = new float[len];
//    color_table_x = &table_color[0];
//    for(int y = 0; y<len; ++y) {
//        *color_table_x++ = exp(-double(y)/sigma_range);
//    }
//    return table_color;
//}

//bilateral filter weights of range
//float * qing_get_range_weighted_table(float sigma_range, int len) {
//    float * table_range, * range_table_x;
//    table_range = new float[len];
//    range_table_x = &table_range[0];
//    for(int y = 0; y<len; ++y) {
//        *range_table_x++ = exp(-double(y*y)/(2*sigma_range*sigma_range));
//    }
//    //cout << endl;
//    return table_range;
//}

////bilateral filter weights of space
//float * qing_get_spatial_weighted_table(float sigma_spatial, int len) {
//    float * table_spatial, * spatial_table_x;
//    table_spatial = new float[len];
//    spatial_table_x = &table_spatial[0];
//    for(int y = 0; y < len; ++y) {
//        *spatial_table_x++ = exp(-double(y*y)/(2*sigma_spatial*sigma_spatial));
//    }
//    // cout << endl;
//    return table_spatial;
//}

//right or wrong
void qing_mcost_to_disp::adaptive_weight_filter(float *out, float *in, uchar * gray_l, uchar * gray_r, int d, int wnd){

    int offset = wnd * 0.5;
    float * temp = new float[m_image_size];
    if(0==temp) {  cerr << "failed to alloc memory.." << endl; return ; }

    float * in_ = in;                  //input
    float * out_ = temp;               //filter result in x-direction

    for(int y = 0; y < m_h; ++y) {
        int idy = y * m_w;
        for(int x = 0; x < m_w; ++x) {
            if(x-d < 0) continue;

            int idx_0 = idy + x;
            int idx_1 = idy + x - d;
            double sum = 0.0, sum_div = 0.0;

            unsigned char * ptr_gray_c_0 = gray_l + idx_0 - 1;
            unsigned char * ptr_gray_c_1 = gray_r + idx_1 - 1;
            unsigned char gray_c_0 = *(++ptr_gray_c_0);         //p
            unsigned char gray_c_1 = *(++ptr_gray_c_1);         //p_d

            cout << "y = " << y << ", x = " << x << endl;
            for(int k = -offset; k <= offset; ++k) {
                if(k+x < 0 || k+x >= m_w || k+x-d<0 || k+x-d >= m_w) continue;

                int idk_0 = idx_0 + k;
                int idk_1 = idx_1 + k;

                unsigned char * ptr_gray_k_0 =  gray_l + idk_0 - 1;    //q
                unsigned char * ptr_gray_k_1 =  gray_r + idk_1 - 1;    //q_d
                int delta_pq = abs(gray_c_0 - *(++ptr_gray_k_0));
                int delta_pq_d = abs(gray_c_1 - *(++ptr_gray_k_1));

                //asw
                //double weight = m_range_table[delta_pq] * m_range_table[delta_pq_d] * m_spatial_table[abs(k)] * m_spatial_table[abs(k)] ;
                //bilateral
                double weight = m_range_table[delta_pq] * m_spatial_table[abs(k)];

                sum_div += weight;
                sum += weight * in_[idk_0];
                cout << "sum = " << sum << ", sum_div = " << sum_div << endl;
            }
            out_[idx_0] = sum / sum_div;
        }
    }
    cout << QING_DEBUG_FLAG_STRING << "\t horizontal filtering done.." << endl;

    //vertical filtering
    in_ = temp;
    out_ = out;
    for(int y = 0; y < m_h; ++y) {
        int idy = y * m_w;
        for(int x = 0; x < m_w; ++x) {
            if(x-d < 0) continue;

            int idx_0 = idy + x;
            int idx_1 = idx_0 - d;
            double sum = 0.0, sum_div = 0.0;

            unsigned char * ptr_gray_c_0 = gray_l + idx_0 - 1;
            unsigned char * ptr_gray_c_1 = gray_r + idx_1 - 1;
            unsigned char gray_c_0 = *(++ptr_gray_c_0);         //p
            unsigned char gray_c_1 = *(++ptr_gray_c_1);         //p_d

            for(int k = -offset; k <= offset; ++k) {
                if(k+y < 0 || k+y >= m_h) continue;
                int idk_0 = idx_0 + k*m_w;
                int idk_1 = idx_1 + k*m_w;

                unsigned char * ptr_gray_k_0 = gray_l + 3 * idk_0 - 1;
                unsigned char * ptr_gray_k_1 = gray_r + 3 * idk_1 - 1;

                int delta_pq = abs(gray_c_0 - *(++ptr_gray_k_0));
                int delta_pq_d = abs(gray_c_1 - *(++ptr_gray_k_1));

                //asw
                //double weight = m_range_table[delta_pq] * m_range_table[delta_pq_d] * m_spatial_table[abs(k)] * m_spatial_table[abs(k)] ;
                //bilateral
                double weight = m_range_table[delta_pq] * m_spatial_table[abs(k)];

                sum_div += weight;
                sum += weight * in_[idk_0];
            }
            out_[idx_0] = sum / sum_div;
        }
    }
    cout << QING_DEBUG_FLAG_STRING << "\t vertical filtering done.." << endl;
}


// sigma_spatial = 0.03; sigma_range = 0.08
void qing_mcost_to_disp::mcost_aggregation(float sigma_range, float sigma_spatial) {

    cout << "start to aggregate matching cost" << endl;

    sigma_range *= QING_FILTER_INTENSITY_RANGE;
    sigma_spatial *= (float)min((float)m_h,(float)m_w);
    //support window
    int wnd = (2*(int)(sigma_spatial+0.5f)) + 1;

    m_range_table = qing_get_range_weighted_table(sigma_range, QING_FILTER_INTENSITY_RANGE);
    m_spatial_table = qing_get_spatial_weighted_table(sigma_spatial, wnd);

    uchar * guidance_l = m_image_l;
    uchar * guidance_r = m_image_r;

    m_filtered_mcost_l = new float[m_total_size];
    if(0==m_filtered_mcost_l) {
        cerr << "failed to initial left filtered matching cost array." << endl;
        exit(-1);
    }
    memcpy(m_filtered_mcost_l, m_mcost_l, sizeof(float)*m_total_size);
    for(int d = 0; d < m_d; ++d) {
        float * out = m_filtered_mcost_l + d * m_image_size;
        adaptive_weight_filter(out, out, guidance_l, guidance_r, d, wnd);  //sigma_spatial == wnd
# if 1
        string out_file = "directional_matching_cost/bf_filtered_mcost_" + qing_int_2_string(d) + ".jpg";
        qing_save_mcost_jpg(out_file, out, m_w, m_h);
# endif
    }

    //filter right-mcost
    //    m_filtered_mcost_r = new float[m_total_size];
    //    if(0==m_filtered_mcost_r) {
    //        cerr << "failed to initial right filtered matching cost array." << endl;
    //        exit(-1);
    //    }
    //    memcpy(m_filtered_mcost_r, m_mcost_r, sizeof(float)*m_total_size);
    //    for(int d = 0; d < m_d; ++d) {
    //        float * out = m_filtered_mcost_r + d * m_image_size;
    //        adaptive_weight_filter(out, out, guidance_l, guidance_r, d, m_h, m_w, wnd, sigma_spatial, sigma_range);
    //    }
}

void qing_mcost_to_disp::directional_adaptive_weight_filter(float *out, float *in, uchar *gray_l, uchar *gray_r, int d, int wnd) {


}


void qing_mcost_to_disp::directional_mcost_aggregation(float sigma_range, float sigma_spatial, const int wnd) {
    cout << "test directional mcost aggregation with zero direction params.." << endl;

    sigma_range *= QING_FILTER_INTENSITY_RANGE;
    sigma_spatial *= min(m_w, m_h);
    //int wnd = 2 * ((int)(sigma_spatial+0.5f)) + 1;
    int offset = wnd * 0.5;
    m_range_table = qing_get_range_weighted_table(sigma_range, QING_FILTER_INTENSITY_RANGE);
    m_spatial_table = qing_get_spatial_weighted_table(sigma_spatial, QING_FILTER_SPATIAL_RANGE);

    cout << "sigma_range = " << sigma_range << endl;
    cout << "sigma_spatial = " << sigma_spatial << endl;
    cout << "wnd = " << wnd << endl;

    uchar * gray_l = m_image_l;
    uchar * gray_r = m_image_r;

    m_filtered_mcost_l = new float[m_total_size];
    if(0==m_filtered_mcost_l) {
        cerr << "failed to new left filtered matching cost..." << endl;
        exit(1);
    }
    std::fill(m_filtered_mcost_l, m_filtered_mcost_l + m_total_size, QING_MAX_MCOST);
    float * min_mcost_x = new float[m_total_size];
    if(0==min_mcost_x) {
        cerr << "failed to new minimum matching cost volume in x-direction..." << endl;
        exit(1);
    }
    std::fill(min_mcost_x, min_mcost_x + m_total_size, QING_MAX_MCOST);

    //generate directions
    int len_wx = 21;
    float wx_step = 1.f/(len_wx/2);
    float * x_directions  =  qing_get_directions(wx_step, len_wx);
# if 1
    for(int len = 0; len < len_wx; ++len) cout << x_directions[len] << ' ';
    cout << endl;
# endif
    int len_wy = 21;
    float wy_step = 1.f/(len_wy/2);
    float * y_directions  =  qing_get_directions(wy_step, len_wy);
# if 1
    for(int len = 0; len < len_wy; ++len) cout << y_directions[len] << ' ';
    cout << endl;
# endif

    double * weights = new double[wnd];       //for each direciton, the weights can be pre-computed
    double sum, sum_div, x_dir, y_dir;
    unsigned char * ptr_gray_c, * ptr_gray_k;
    unsigned char gray_c;
    int idy, idx, idk, delta_pq;

    QingTimer timer;

    //x-direction filtering
    cout << "horizontal filtering start.." << endl;

    for(int d = 0; d < m_d; ++d) {
        float * out = min_mcost_x + d * m_image_size;
        cout << "d = " << d ; timer.restart();

        for(int y = 0; y < m_h; ++y) {
            idy = y * m_w;
            for(int x = 0; x < m_w; ++x) {
                idx = idy + x;

                ptr_gray_c = gray_l + idx - 1;
                gray_c = *(++ptr_gray_c);

                std::fill(weights, weights + wnd, 0.0);
                for(int k = -offset, tk = 0; k <= offset; ++tk, ++k) {
                    if(x+k < 0 || x+k >= m_w) continue;
                    idk = idx + k;

                    ptr_gray_k = gray_l + idk - 1;
                    delta_pq = abs(gray_c - *(++ptr_gray_k));
                    weights[tk] = m_range_table[delta_pq] * m_spatial_table[abs(k)];
                }

                //compute aggregate matching cost in this horizontal support window of each directions
                //select the minimum , stored in out[idx], i.e, min_mcost_x[d*m_image_size+idx]
                for(int wx = 0; wx < len_wx; ++wx) {
                    x_dir = x_directions[wx];
                    sum = 0.0;
                    sum_div = 0.0;

                    for(int k = -offset, tk = 0; k<=offset; ++tk,++k) {
                        if(x+k < 0 || x+k >= m_w) continue;
                        idk = idx + k;

                        int k_d = d + x_dir * k;
//                        if(k_d < 0 || k_d >= m_d) continue;
                        k_d = max(0,k_d);
                        k_d = min(k_d, m_d-1);

                        sum += weights[tk] * (*(m_mcost_l + k_d * m_image_size + idk));
                        sum_div += weights[tk];
                        //                        if(x == 30)
                        //                            cout << "x = " << x << " sum = " << sum << ", sum_div = " << sum_div << endl;
                    }
                    //if(x==30)  exit(1);
                    if(sum_div >= 0.000001) sum/=sum_div;
                    if(sum < out[idx]) { out[idx] = sum; }
                }
            }
        }
        cout << "\t duration = " << timer.duration() * 1000 << " ms." << endl;

# if 0
        qing_create_dir("aggr_matching_cost");
        string out_file = "aggr_matching_cost/x_mcost_" + qing_int_2_string(d) + ".jpg";
        qing_save_mcost_jpg(out_file, out, m_w, m_h);
        string out_txt_file = "aggr_matching_cost/x_mcost_" + qing_int_2_string(d) + ".txt";
        qing_save_mcost_txt(out_txt_file, out, m_image_size, m_w);
# endif
    }

    cout << "horizontal filtering end..." << endl;

    //y-directoin filtering
    cout << "vertical filtering start..." << endl;
    for(int d = 0; d < m_d; ++d) {
        float * out = m_filtered_mcost_l + d * m_image_size;
        cout << "d = " << d ; timer.restart();

        for(int y = 0; y < m_h; ++y) {
            idy = y * m_w;
            for(int x = 0; x < m_w; ++x) {
                idx = idy + x;

                ptr_gray_c = gray_l + idx - 1;
                gray_c = *(++ptr_gray_c);

                std::fill(weights, weights+wnd, 0.0);
                for(int k = -offset, tk = 0; k <= offset; ++tk, ++k) {
                    if(y+k < 0 || y+k >= m_h) continue;

                    idk = idx + k * m_w;
                    ptr_gray_k = gray_l + idk - 1;
                    delta_pq = abs(gray_c - *(++ptr_gray_k));
                    weights[tk] = m_range_table[delta_pq] * m_spatial_table[abs(k)];
                }

                //compute aggregate matching cost in this vertical support window of each directions
                //select the minumum, stored in out[idx], i.e, m_filtered_mcost_l[d*m_image_size+idx]

                for(int wy = 0; wy < len_wy; ++wy) {
                    y_dir = y_directions[wy];
                    sum = 0.0;
                    sum_div = 0.0;
                    for(int k = -offset, tk = 0; k<=offset; ++tk, ++k) {
                        if(y+k < 0 || y+k >= m_w) continue;
                        idk = idx + k * m_w;
                        int k_d = d + y_dir * k;
                       // if(k_d < 0 || k_d >= m_d) continue;
                        k_d = max(0,k_d);
                        k_d = min(k_d, m_d-1);

                        sum += weights[tk] * (*(min_mcost_x + k_d * m_image_size + idk));
                        sum_div += weights[tk];
                    }

                    if(sum_div >= 0.0000001) sum/=sum_div;
                    if(sum < out[idx]) {out[idx] = sum; }
                }
            }
        }
          cout << "\t duration = " << timer.duration() * 1000 << " ms." << endl;
# if 0
        qing_create_dir("aggr_matching_cost");
        string out_file = "aggr_matching_cost/y_mcost_" + qing_int_2_string(d) + ".jpg";
        qing_save_mcost_jpg(out_file, out, m_w, m_h);
        string out_txt_file = "aggr_matching_cost/y_mcost_" + qing_int_2_string(d) + ".txt";
        qing_save_mcost_txt(out_txt_file, out, m_image_size, m_w);
# endif
    }

    cout << "vertial filtering end...." << endl;
}


void qing_mcost_to_disp::mcost_to_disp(const int scale, const string savename) {
    memset(m_disp_l, 0, sizeof(unsigned char) * m_image_size);
    memset(m_disp_r, 0, sizeof(unsigned char) * m_image_size);
    for(int y = 0, idx = 0; y < m_h; ++y) {
        for(int x = 0; x < m_w; ++x) {
            float min_mcost = 2 * QING_MAX_MCOST;
            int min_d = 0;
            for(int d = 1; d < m_d; ++d) {
                float mcost = m_filtered_mcost_l[d*m_image_size + idx];
                if(mcost < min_mcost) {
                    min_mcost = mcost;
                    min_d = d;

                }
            }
            m_disp_l[idx] = (unsigned char)(min_d * scale);
            idx ++;
        }
    }
    //    for(int y = 0, idx = 0; y < m_h; ++y) {
    //        for(int x = 0; x < m_w; ++x) {
    //            float min_mcost = 10000.f;
    //            float min_d = 0;
    //            for(int d = 1; d < m_d; ++d) {
    //                float mcost = m_filtered_mcost_r[d*m_image_size + idx];
    //                if(mcost < min_mcost) {
    //                    min_mcost = mcost;
    //                    min_d = d;
    //                }
    //            }
    //            m_disp_r[idx] = min_d * scale;
    //        }
    //    }

# if 1
    Mat left_disp(m_h, m_w, CV_8UC1);
    memcpy(left_disp.data, m_disp_l, sizeof(unsigned char)*m_image_size);
    Mat right_disp(m_h, m_w, CV_8UC1);
    memcpy(right_disp.data, m_disp_r, sizeof(unsigned char)*m_image_size);
    cout << "save ./left_" << savename << endl; imwrite("./left_" + savename,  left_disp);
//    cout << "save ./disp_" << savename << endl; imwrite("./disp_" + savename,  left_disp);
//    cout << "save ./right_" << savename << endl; imwrite("./right_" + savename, right_disp);
# endif
}

