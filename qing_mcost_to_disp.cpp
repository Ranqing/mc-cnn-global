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

void qing_save_mcost_vol(const string filename, float * cost_vol, int m_d, int m_h, int m_w) {
    fstream fout(filename.c_str(), ios::out);
    if(fout.is_open() == false) {
        cerr << "failed to open " << filename << endl;
        return ;
    }
    cout << "saving " << filename << endl;
    int m_image_size = m_h * m_w;
    for(int d = 0; d < m_d; ++d) {
        float * mcost = cost_vol + d * m_image_size;
        for(int y = 0, idx = 0; y < m_h; ++y) {
            for(int x = 0; x < m_w; ++x) {
                fout << mcost[idx++] << ' ';
            }
            fout << endl;
        }
        cout << d << ' ';
    }
    cout << endl;    fout.close();
}

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
void qing_mcost_to_disp::mcost_aggregation(const int wnd) {

    cout << "start to aggregate matching cost" << endl;
    cout << "filtering wnd = " << wnd << endl;

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


void qing_mcost_to_disp::directional_mcost_aggregation(const int wnd) {
    // cout << "test directional mcost aggregation with zero direction params.." << endl;
    cout << "directional mcost aggregation with direction params.." << endl;
    cout << "wnd = " << wnd << endl;

    uchar * gray_l = m_image_l;
    uchar * gray_r = m_image_r;

    int len = 21;
    float step = 1.0f/(len*0.5f);
    float * directions = qing_get_directions(step, len);

    m_filtered_mcost_l = new float[m_total_size];
    if(0==m_filtered_mcost_l) {
        cerr << "failed to new left filtered matching cost..." << endl;
        exit(1);
    }
    std::fill(m_filtered_mcost_l, m_filtered_mcost_l + m_total_size, QING_MAX_MCOST);
    float * min_mcost_x = new float[m_total_size];
    if(0==min_mcost_x) {
        cerr << "failed to new minimum x matcihng cost..." << endl;
        exit(1);
    }
    std::fill(min_mcost_x, min_mcost_x + m_total_size, QING_MAX_MCOST);
    qing_directional_bf_mcost_aggregation(m_filtered_mcost_l, m_mcost_l, min_mcost_x, gray_l, m_w, m_h, m_d, wnd, m_range_table, m_spatial_table, directions, len);

# if STEREO_RIGHT
    m_filtered_mcost_r = new float[m_total_size];
    if(0==m_filtered_mcost_r) {
        cerr << "failed to new right filtered matching cost..." << endl;
        exit(1);
    }
    std::fill(m_filtered_mcost_r, m_filtered_mcost_r + m_total_size, QING_MAX_MCOST);
    std::fill(min_mcost_x, min_mcost_x + m_total_size, QING_MAX_MCOST);
    qing_directional_bf_mcost_aggregation(m_filtered_mcost_r, m_mcost_r, min_mcost_x, gray_r, m_w, m_h, m_d, wnd, m_range_table, m_spatial_table, directions, len);
# endif
}

void qing_mcost_to_disp::save_filtered_mcost(const string folder) {
    string filename = folder + "/left.txt" ;
    qing_save_mcost_vol(filename, m_filtered_mcost_l, m_d, m_h, m_w);
    filename = folder + "/left.bin";
    qing_write_bin(filename, m_filtered_mcost_l, m_total_size);

# if STEREO_RIGHT
    filename = folder + "/right.txt";
    qing_save_mcost_vol(filename, m_filtered_mcost_r, m_d, m_h, m_w);
    filename = folder + "/right.bin";
    qing_write_bin(filename, m_filtered_mcost_r, m_total_size);
# endif
}


void qing_mcost_to_disp::mcost_to_disp(const int scale, const string savename) {
    memset(m_disp_l, 0, sizeof(unsigned char) * m_image_size);
    for(int y = 0, idx = 0; y < m_h; ++y) {
        for(int x = 0; x < m_w; ++x) {
            float min_mcost = 2 * QING_MAX_MCOST;
            int min_d = 0;
            for(int d = 1; d < m_d; ++d) {
                //   float mcost = m_filtered_mcost_l[d*m_image_size + idx];
                float mcost = m_mcost_l[d*m_image_size + idx];
                if(mcost < min_mcost) {
                    min_mcost = mcost;
                    min_d = d;

                }
            }
            m_disp_l[idx] = (unsigned char)(min_d * scale);
            idx ++;
        }
    }
    Mat left_disp(m_h, m_w, CV_8UC1);
    memcpy(left_disp.data, m_disp_l, sizeof(unsigned char)*m_image_size);
    cout << "save ./left_" << savename << endl; imwrite("./left_" + savename,  left_disp);

# if STEREO_RIGHT
    memset(m_disp_r, 0, sizeof(unsigned char) * m_image_size);

    for(int y = 0, idx = 0; y < m_h; ++y) {
        for(int x = 0; x < m_w; ++x) {
            float min_mcost = 10000.f;
            float min_d = 0;
            for(int d = 1; d < m_d; ++d) {
                //                float mcost = m_filtered_mcost_r[d*m_image_size + idx];
                float mcost = m_mcost_r[d*m_image_size + idx];
                if(mcost < min_mcost) {
                    min_mcost = mcost;
                    min_d = d;
                }
            }
            m_disp_r[idx] = min_d * scale;
        }
    }
    Mat right_disp(m_h, m_w, CV_8UC1);
    memcpy(right_disp.data, m_disp_r, sizeof(unsigned char)*m_image_size);
    cout << "save ./right_" << savename << endl; imwrite("./right_" + savename, right_disp);
# endif

}

