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

//bilateral filter weights
float * qing_get_range_weighted_table(float sigma_range, int len) {
    float * table_range, * range_table_x;
    table_range = new float[len];
    range_table_x = &table_range[0];
    for(int y = 0; y<len; ++y) {
        *range_table_x++ = exp(-double(y*y)/(2*sigma_range*sigma_range));
    }
    //cout << endl;
    return table_range;
}

//can be seperated into two-dimensional
float * qing_get_spatial_weighted_table(float sigma_spatial, int len) {
    float * table_spatial, * spatial_table_x;
    table_spatial = new float[len];
    spatial_table_x = &table_spatial[0];
    for(int y = 0; y < len; ++y) {
        *spatial_table_x++ = exp(-double(y*y)/(2*sigma_spatial*sigma_spatial));
    }
   // cout << endl;
    return table_spatial;
}

void qing_mcost_to_disp::adaptive_weight_filter(float *out, float *in, uchar *color_0, uchar *color_1, int d, int h, int w, float sigma_spatial, float sigma_range) {

    float * temp = new float[h*w];
    if(0==temp) {
        cerr << "failed to alloc memory.." << endl;
        return ;
    }
    float * in_ = in;                  //input
    float * out_ = temp;               //filter result in x-direction

    int wnd_offset = (int)(sigma_spatial * 0.5 + 0.5f);

    //gudiance filter weight can be computed first : no since the memory cost to sigma_range * sigma_range * w * h
    //brute force: result of baseline

    for(int y = 0; y < h; ++y) {
        int idy = y * w;
        for(int x = 0; x < d; ++x) {
            int idx = idy + x;
            out_[idx] = in_[idx];
        }
        for(int x = d; x < w; ++x) {
            int idx_0 = idy + x;                 //center pixel
            int idx_1 = idx_0 - d;               //corresponding pixel

            float res = 0.0f;

            //can be move outside
//            uchar * ptr_c_bgr_0 = bgr_0+m_channels*idx_0-1;      //color for center pixel
//            uchar * ptr_c_bgr_1 = bgr_1+m_channels*idx_1-1;      //color for corresponding pixel
//            uchar c_b_0, c_g_0, c_r_0, c_b_1, c_g_1, c_r_1;
//            c_b_0 = *(++ptr_c_bgr_0); c_g_0 = *(++ptr_c_bgr_0); c_r_0 = *(++ptr_c_bgr_0);
//            c_b_1 = *(++ptr_c_bgr_1); c_g_1 = *(++ptr_c_bgr_1); c_r_1 = *(++ptr_c_bgr_1);
            uchar * ptr_c_color_0 = color_0 + m_channels*idx_0-1;
            uchar * ptr_c_color_1 = color_1 + m_channels*idx_1-1;
            uchar c_color_0 = *(++ptr_c_color_0);
            uchar c_color_1 = *(++ptr_c_color_1);

            cout << '.';

            for(int j = -wnd_offset; j <= wnd_offset; ++j) {
                if(0>j+y||h<=j+y) continue;
                for(int i = -wnd_offset; i <= wnd_offset; ++i) {
                    if(0>i+x||w<=i+x) continue;
                    if(0>i+x-d) continue;

                    int k_idx_0 = idx_0 + j*w+i;
                    int k_idx_1 = k_idx_0 - d;

//                    uchar * ptr_k_bgr_0 = bgr_0 + 3 * k_idx_0 - 1;
//                    uchar * ptr_k_bgr_1 = bgr_1 + 3 * k_idx_1 - 1;

//                    int delta_0 = (int)(abs(*(++ptr_k_bgr_0) - c_b_0));
//                    int delta_1 = (int)(abs(*(++ptr_k_bgr_0) - c_g_0));
//                    int delta_2 = (int)(abs(*(++ptr_k_bgr_0) - c_r_0));
//                    int delta_3 = (int)(abs(*(++ptr_k_bgr_1) - c_b_1));
//                    int delta_4 = (int)(abs(*(++ptr_k_bgr_1) - c_g_1));
//                    int delta_5 = (int)(abs(*(++ptr_k_bgr_1) - c_r_1));

//                    int c_diff_0 = (int)(sqrt((delta_0 * delta_0 + delta_1 * delta_1 + delta_2 * delta_2)) + .5f);
//                    int c_diff_1 = (int)(sqrt((delta_3 * delta_3 + delta_4 * delta_4 + delta_5 * delta_5)) + .5f);

                    uchar * ptr_k_color_0 = color_0 + k_idx_0 -1;
                    uchar * ptr_k_color_1 = color_1 + k_idx_1 -1;
                    int delta0 = (int)(abs(*(++ptr_k_color_0) - c_color_0));
                    int delta1 = (int)(abs(*(++ptr_k_color_1) - c_color_1));
                    float spatial_w = m_spatial_table[abs(i)] * m_spatial_table[abs(j)];
                    double weight = m_range_table[delta0] * m_range_table[delta1] * spatial_w * spatial_w;

                 //   cout << "weight = " << weight << endl;
                    res += weight * in_[k_idx_0];
                }
            }
            out_[idx_0] = res;
        }
    }
    cout << endl;
    memcpy(out, out_, sizeof(float)*m_image_size);
   // getchar();
    cout << "COOL" << endl;
}


// sigmaSpatial defaults to min( width, height )  * 0.03
// sigmaRange defaults to ( edgeMax - edgeMin ) * 0.08
// or sigma_spatial = 0.03; sigma_range = 0.08
void qing_mcost_to_disp::mcost_aggregation(float sigma_range, float sigma_spatial) {
    cout << "start to aggregate matching cost" << endl;

    sigma_spatial *= (float)min((float)m_h,(float)m_w);
    sigma_range *= QING_FILTER_INTENSITY_RANGE;
    cout << sigma_spatial << ", " << sigma_range << endl;
    m_range_table = qing_get_range_weighted_table(sigma_range, QING_FILTER_INTENSITY_RANGE);
    m_spatial_table = qing_get_spatial_weighted_table(sigma_spatial, (int)(sigma_spatial+0.5f));
    cout << "check range table..." << endl;
    for(int i = 0; i < QING_FILTER_INTENSITY_RANGE; ++i) {
        cout << m_range_table[i] << ' ' ;
    }
    cout << endl;
    cout << "check spatial table..." << endl;
    for(int i = 0; i < QING_FILTER_INTENSITY_RANGE; ++i) {
        cout << m_spatial_table[i] << ' ' ;
    }
    cout << endl;
    cout << "initialization of range_table and spatial_table done.." << endl;

    uchar * guidance_l = m_image_l;
    uchar * guidance_r = m_image_r;
    float * ones = new float[m_image_size];

    memcpy(m_filtered_mcost_l, m_mcost_l, sizeof(float)*m_total_size);
    for(int d = 0; d < m_d; ++d) {
        memset(ones, 1, sizeof(float)*m_image_size);

        float * out = m_filtered_mcost_l + d * m_image_size;

        adaptive_weight_filter(out,  out,  guidance_l, guidance_r, d, m_h, m_w, sigma_spatial, sigma_range);  //sigma_spatial == wnd
        adaptive_weight_filter(ones, ones, guidance_l, guidance_r, d, m_h, m_w, sigma_spatial, sigma_range);

        for(int i = 0; i < m_image_size; ++i) {
            if(ones[i] > 0.00001f) out[i] = out[i]/ones[i];
        }
    }

    //filter right-mcost
    memcpy(m_filtered_mcost_r, m_mcost_r, sizeof(float)*m_total_size);
    for(int d = 0; d < m_d; ++d) {
        memset(ones, 1, sizeof(float)*m_image_size);

        float * out = m_filtered_mcost_r + d * m_image_size;

        adaptive_weight_filter(out, out, guidance_l, guidance_r, d, m_h, m_w, sigma_spatial, sigma_range);
        adaptive_weight_filter(out, out, guidance_r, guidance_r, d, m_h, m_w, sigma_spatial, sigma_range);
        for(int i = 0; i < m_image_size; ++i) {
            if(ones[i] > 0.00001f) out[i] = out[i]/ones[i];
        }
    }
}

void qing_mcost_to_disp::mcost_to_disp() {

    for(int j = 0; j < m_h; ++j) {
        for(int i = 0; i < m_w; ++i) {
            int idx = j * m_w + i;
            int max_d_l = 0;
            float min_mcost_l = 1;
            int max_d_r = 0;
            float min_mcost_r = 1;

            for(int d = 0; d < m_d; ++d) {
                int d_idx = d * m_image_size + idx;
                float d_mcost_l = m_filtered_mcost_l[d_idx];
                float d_mcost_r = m_filtered_mcost_r[d_idx];

                if(false == qing_is_nan(d_mcost_l)) {
                    if(d_mcost_l < min_mcost_l) {
                        min_mcost_l = d_mcost_l;
                        max_d_l = d;
                    }
                }
                if(false == qing_is_nan(d_mcost_r)) {
                    if(d_mcost_r < min_mcost_r) {
                        min_mcost_r = d_mcost_r;
                        max_d_r = d;
                    }
                }
            }
            m_disp_l[idx] = max_d_l;
            m_disp_r[idx] = max_d_r;
        }
    }

# if 1
    Mat left_disp(m_h, m_w, CV_8UC1);
    uchar * data = left_disp.ptr<uchar>(0);
    for(int j = 0, size = 0; j < m_h; ++j) {
        for(int i = 0; i < m_w; ++i) {
            data[size] = m_disp_l[size];
            size ++;
        }
    }
    Mat right_disp(m_h, m_w, CV_8UC1);
    data = right_disp.ptr<uchar>(0);
    for(int j = 0, size = 0; j < m_h; ++j) {
        for(int i = 0; i < m_w; ++i) {
            data[size] = m_disp_r[size];
            size ++;
        }
    }
//    imshow("left_disp", left_disp );
//    imshow("right_disp", right_disp );
//    waitKey(0);
//    destroyAllWindows();
    imwrite("./left_filtered.png", left_disp);
    imwrite("./disp_filtered.png", left_disp);
    imwrite("./right_filtered.png", right_disp);
# endif
}

