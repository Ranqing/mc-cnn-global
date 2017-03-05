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

            for(int k = -offset; k <= offset; ++k) {
                if(k+x < 0 || k+x >= m_w || k+x-d<0 || k+x-d >= m_w) continue;

                int idk_0 = idx_0 + k;
                int idk_1 = idx_1 + k;

                unsigned char * ptr_gray_k_0 =  gray_l + idk_0 - 1;    //q
                unsigned char * ptr_gray_k_1 =  gray_r + idk_1 - 1;    //q_d
                int delta_pq = abs(gray_c_0 - *(++ptr_gray_k_0));
                int delta_pq_d = abs(gray_c_1 - *(++ptr_gray_k_1));

                //asw
                double weight = m_range_table[delta_pq] * m_range_table[delta_pq_d] * m_spatial_table[abs(k)] * m_spatial_table[abs(k)] ;
                //bilateral
                //double weight = m_range_table[delta_pq] * m_spatial_table[abs(k)];

                sum_div += weight;
                sum += weight * in_[idk_0];
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

                double weight = m_range_table[delta_pq] * m_range_table[delta_pq_d] * m_spatial_table[abs(k)] * m_spatial_table[abs(k)] ;
                //bilateral
                //double weight = m_range_table[delta_pq] * m_spatial_table[abs(k)];

                sum_div += weight;
                sum += weight * in_[idk_0];
            }
            out_[idx_0] = sum / sum_div;
        }
    }
    cout << "\t vertical filtering done.." << endl;
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
# if 0
        string out_file = "matching_cost/filtered_mcost_" + qing_int_2_string(d) + ".jpg";
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

void qing_mcost_to_disp::directional_mcost_aggregation(float sigma_range, float sigma_spatial) {

    cout << "start to aggregate matching cost with direction params" << endl;
    sigma_range *= QING_FILTER_INTENSITY_RANGE;
    sigma_spatial *= (float)min((float)m_h,(float)m_w);
    //support window
    int wnd = (2*(int)(sigma_spatial+0.5f)) + 1;
    int offset = wnd * 0.5;

    m_range_table = qing_get_range_weighted_table(sigma_range, QING_FILTER_INTENSITY_RANGE);
    m_spatial_table = qing_get_spatial_weighted_table(sigma_spatial, wnd);

    uchar * gray_l = m_image_l;
    uchar * gray_r = m_image_r;

    //final filtered results
    cout << m_filtered_mcost_l << endl;
    m_filtered_mcost_l = new float[m_total_size];
    if(0==m_filtered_mcost_l) {
        cerr << "failed to initial left filtered matching cost array." << endl;
        exit(-1);
    }
    //   memset(m_filtered_mcost_l, QING_MAX_MCOST, sizeof(float)*m_total_size); set value each byte
    std::fill(m_filtered_mcost_l, m_filtered_mcost_l + m_total_size, QING_MAX_MCOST);


    //minimum matching cost in x-direction
    float * min_mcost_x = new float[m_total_size];
    if(0==min_mcost_x) {
        cerr << "failed to initial intermediate matching cost array in x-direction" << endl;
        exit(-1);
    }
    //    memset(min_mcost_x, QING_MAX_MCOST, sizeof(float)*m_total_size);
    std::fill(min_mcost_x, min_mcost_x + m_total_size, QING_MAX_MCOST);

    double x_directions[5] = {-0.02, -0.01, 0, 0.01, 0.02};
    double y_directions[5] = {-0.02, -0.01, 0, 0.01, 0.02};

    //x-direction
    //in: m_mcost_l
    //out:min_mcost_x
    cout << QING_DEBUG_FLAG_STRING << "\t horizontal filtering start.." << endl;
    double * weights = new double[wnd+4];

    for(int d = 0; d < m_d; ++d) {
        float * out = min_mcost_x + d * m_image_size;
        cout << "d = " << d << endl;

        for(int y = 0; y < m_h; ++y) {
            int idy = y * m_w;
            for(int x = 0; x < m_w; ++x) {
                int idx = idy + x;

                //compute aggregate x-directional minimum matching cost: min_mcost_x[d*m_image_size + idx];

                unsigned char * ptr_gray_c = gray_l + idx - 1;
                unsigned char gray_c = *(++ptr_gray_c);

                //pre-compute weights for the neighbors
                memset(weights, 0, sizeof(double)*(wnd+4));
                for(int k = -offset, tk = 0; k <= offset; ++tk, ++k) {
                    if(x+k < 0 || x+k >= m_w) continue;
                    int idk = idx + k;
                    unsigned char * ptr_gray_k = gray_l + idk - 1;
                    int delta_pq = abs(gray_c - *(++ptr_gray_k));
                    weights[tk] = m_range_table[delta_pq] * m_spatial_table[abs(k)];
                }

                //aggregate initial matching cost in x-direction with a direction param
                //select the minimum x-aggregated matching cost, saving in out[idx]
                for(int wx = 0; wx < 5; ++wx) {
                    float x_dir = x_directions[wx];
                    float sum = 0.f;
                    float sum_div = 0.f;
                    for(int k = -offset, tk = 0; k <=offset; ++tk, ++k) {
                        int k_d = d + x_dir * k;
                        int idk = idx + k;
                        if(k_d < 0 || k_d >= m_d) continue;

                        sum += weights[tk] * (*(m_mcost_l + k_d * m_image_size + idk));
                        sum_div += weights[tk];
                    }
                    if(sum_div > 0.0000001)  sum /= sum_div;
                    if(out[idx] > sum) out[idx] = sum;
                }
            }
        }
    }
    cout << QING_DEBUG_FLAG_STRING << "\t horizontal filtering done.." << endl;


    //y-direction
    //in: min_mcost_x
    //out: m_filtered_mcost_x
    cout << QING_DEBUG_FLAG_STRING << "\t vertical filtering start.." << endl;
    for(int d = 0; d < m_d; ++d) {
        float * out = m_filtered_mcost_l + d * m_image_size;
        cout << "d = " << d << endl;

        for(int y = 0; y < m_h; ++y) {
            int idy = y * m_w;
            for(int x = 0; x < m_w; ++x) {
                int idx = idy + x;

                //compute aggregate y-directional minimum matching cost: m_filtered_mcost[d*m_image_size+idx]

                unsigned char * ptr_gray_c = gray_l + idx - 1;
                unsigned char gray_c = *(++ptr_gray_c);

                //pre_compute weights for neigbhors
                memset(weights, 0, sizeof(double)*(wnd+4));
                for(int k = -offset, tk = 0; k <= offset; ++tk, ++k) {
                    if(y+k < 0 || y+k>=m_h) continue;
                    int idk = idx + k * m_w;

                    unsigned char * ptr_gray_k = gray_l + idk - 1;
                    int delta_pq = abs(gray_c - *(++ptr_gray_k));
                    weights[tk] = m_range_table[delta_pq] * m_spatial_table[abs(k)];
                }

                //aggregate initial matching cost in y-direction with a direction param
                //select the minimum y-aggegate matching cost, saving in out[idx]
                for(int wy = 0; wy < 5; ++wy) {
                    float y_dir = y_directions[wy];
                    float sum = 0.0f;
                    float sum_div = 0.0f;
                    for(int k = -offset, tk = 0; k <= offset; ++tk, ++k) {
                        if(y+k < 0 || y+k>=m_h) continue;
                        int idk = idx + k * m_w;
                        int k_d = d + y_dir * k;
                        if(k_d < 0 || k_d >= m_d) continue;

                        sum += weights[tk] * (*(min_mcost_x + k_d * m_image_size + idk));
                        sum_div += weights[tk];
                    }

                    if(sum_div > 0.000001) sum /= sum_div;
                    if(out[idx] > sum) out[idx] = sum;
                }
            }
        }
# if 1
        qing_create_dir("./directional_matching_cost/");
        string out_file = "./directional_matching_cost/directional_mcost_" + qing_int_2_string(d) + ".jpg";
        qing_save_mcost_jpg(out_file, out, m_w, m_h);
# endif
    }
    cout << QING_DEBUG_FLAG_STRING << "\t vertical filtering done.." << endl;

}


void qing_mcost_to_disp::mcost_to_disp(const int scale, const string savename) {
    memset(m_disp_l, 0, sizeof(unsigned char) * m_image_size);
    memset(m_disp_r, 0, sizeof(unsigned char) * m_image_size);
    for(int y = 0, idx = 0; y < m_h; ++y) {
        for(int x = 0; x < m_w; ++x) {
            float min_mcost = 2 * QING_MAX_MCOST;
            float min_d = 0;
            for(int d = 1; d < m_d; ++d) {
                float mcost = m_filtered_mcost_l[d*m_image_size + idx];
                if(mcost < min_mcost) {
                    min_mcost = mcost;
                    min_d = d;
                }
            }
            m_disp_l[idx] = min_d * scale;
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
    imwrite("./left_" + savename,  left_disp);
    imwrite("./disp_" + savename,  left_disp);
    imwrite("./right_" + savename, right_disp);
# endif
}

