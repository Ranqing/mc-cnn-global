#include "qing_mcost_to_disp.h"
#include "../../../Qing/qing_scanline_optimizer.h"

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

    uchar * gray_l = m_image_l;
    uchar * gray_r = m_image_r;

    m_filtered_mcost_l = new float[m_total_size];
    if(0==m_filtered_mcost_l) {
        cerr << "failed to initial left filtered matching cost vol..." << endl;
        exit(-1);
    }
    cout << "left" << endl;
    qing_bf_mcost_aggregation(m_filtered_mcost_l, m_mcost_l, gray_l, m_w, m_h, m_d, wnd, m_range_table, m_spatial_table);

# if STEREO_RIGHT
    m_filtered_mcost_r = new float[m_total_size];
    if(0==m_filtered_mcost_r) {
        cerr << "failed to initial right filtered matching cost vol..." << endl;
        exit(-1);
    }
    cout << "right" << endl;
    qing_bf_mcost_aggregation(m_filtered_mcost_r, m_mcost_r, gray_r, m_w, m_h, m_d, wnd, m_range_table, m_spatial_table);
# endif
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
 //   qing_directional_aw_mcost_aggregation_l(m_filtered_mcost_l, m_mcost_l, min_mcost_x, gray_l, gray_r, m_w, m_h, m_d, wnd, m_range_table, m_spatial_table, directions, len);

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

void qing_mcost_to_disp::mcost_to_disp(const int scale, const string savename) {
    memset(m_disp_l, 0, sizeof(unsigned char) * m_image_size);
    for(int y = 0, idx = 0; y < m_h; ++y) {
        for(int x = 0; x < m_w; ++x) {
            float min_mcost = 2 * QING_MAX_MCOST;
            int min_d = 0;
            for(int d = 1; d < m_d; ++d) {
                float mcost = m_filtered_mcost_l[d*m_image_size + idx];
                //  float mcost = m_mcost_l[d*m_image_size + idx];
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
                float mcost = m_filtered_mcost_r[d*m_image_size + idx];
                // float mcost = m_mcost_r[d*m_image_size + idx];
                if(mcost < min_mcost) {
                    min_mcost = mcost;
                    min_d = d;
                }
            }
            m_disp_r[idx] = min_d * scale;
            idx++;
        }
    }
    Mat right_disp(m_h, m_w, CV_8UC1);
    memcpy(right_disp.data, m_disp_r, sizeof(unsigned char)*m_image_size);
    cout << "save ./right_" << savename << endl; imwrite("./right_" + savename, right_disp);
# endif
}

void qing_mcost_to_disp::semi_global() {

}

// update cost with scanline optimization method:
// Cr(p,d) = C1(p,d) + min( Cr(p-r,d), Cr(p-r, d+1)+P1,  Cr(p-r, d-1)+P1, min_k(Cr(p-r,k)+P2) ) - min_k(Cr(p-r,k))
void qing_mcost_to_disp::scanline_optim(float *so_mcost_vol, float *in_mcost_vol, int x1, int x2, int y1, int y2) {

    int idx1 = y1 * m_w + x1;
    int idx2 = y2 * m_w + x2;      //previous pixel in direction

//  cout << x1 << ' ' << x2 << ' ' << y1 << ' ' << y2 <<  endl;
//    bool flag = false;
//    if(x1 == 116 && x2 == 116 && y1 == 185 && y2 == 184) {
//        flag = true;
//        cout << "HEHE" << endl;
//    }
//    if(x1 == 193 && x2 == 193 && y1 == 56 && y2 == 55) {
//        flag = true;
//        cout << "HEHE" << endl;
//    }
//    if(x1 == 0 && x2 == 0 && y1 == 183 && y2 == 184) {
//        flag = true;
//        cout << "HEHE" << endl;
//    }

    //i ~ min[L_r(p_r, i)]
    float * ptr_mcost = so_mcost_vol;
    float min_so_mcost = ptr_mcost[idx2];
    for(int d = 1; d < m_d; ++d) {
        ptr_mcost += m_image_size;
        if(ptr_mcost[idx2] < min_so_mcost)
            min_so_mcost = ptr_mcost[idx2];
    }
    //if(flag) cout << min_so_mcost << endl;

    float min_k_cr = min_so_mcost;
    float * ptr_in_mcost_vol = in_mcost_vol - m_image_size;
    float * ptr_so_mcost_vol = so_mcost_vol - m_image_size;
    float * ptr_add_so_mcost_vol, * ptr_sub_so_mcost_vol;
    float so_mcost, tmp_mcost;
    float P1, P2;
    for(int d = 0; d < m_d; ++d) {
        ptr_in_mcost_vol +=  m_image_size;
        ptr_so_mcost_vol += m_image_size;

        so_mcost = *(ptr_in_mcost_vol + idx1) - min_k_cr;

        scanline_optim_params(x1, x2, y1, y2, d, P1, P2);
 //       if(flag) cout  << "d = " << d << "\tP1 = " << P1 << "\tP2 = " << P2 << endl;

        min_so_mcost = min_k_cr + P2;
        if(min_so_mcost > ptr_so_mcost_vol[idx2])
            min_so_mcost = ptr_so_mcost_vol[idx2];

        if( d != 0) {
            ptr_sub_so_mcost_vol = ptr_so_mcost_vol - m_image_size;
            if( min_so_mcost > (tmp_mcost = ptr_sub_so_mcost_vol[idx2] + P1) )
                min_so_mcost = tmp_mcost;
        }
        if( d != m_d-1) {
            ptr_add_so_mcost_vol = ptr_so_mcost_vol + m_image_size;
            if( min_so_mcost > (tmp_mcost = ptr_add_so_mcost_vol[idx2] + P1) )
                min_so_mcost = tmp_mcost;
        }

        ptr_so_mcost_vol[idx1] = so_mcost + min_so_mcost;
    }
}

//pi1 = 1.0;
//pi2 = 3.0;
//tauSO = 15;
void qing_mcost_to_disp::scanline_optim_params(int x1, int x2, int y1, int y2, int d, float &P1, float &P2) {
    int idx1 = y1 * m_w + x1;
    int idx2 = y2 * m_w + x2;
    int diff1 = abs(m_image_l[idx2] - m_image_l[idx1]);
    int diff2 = QING_SO_TAU + 1;
    if( 0<=x1-d && m_w>x1-d && 0<=x2-d && m_w>x2-d) {
        diff2 = abs(m_image_r[idx2-d] - m_image_r[idx1-d]);
    }
    if(diff1 < QING_SO_TAU) {
        if(diff2 < QING_SO_TAU) {
            P1 = QING_SO_P1;
            P2 = QING_SO_P2;
        }
        else {
            P1 = float(QING_SO_P1/4.0);
            P2 = float(QING_SO_P2/4.0);
        }
    }
    else{
        if(diff2 < QING_SO_TAU) {
            P1 = float(QING_SO_P1/4.0);
            P2 = float(QING_SO_P2/4.0);
        }
        else {
            P1 = float(QING_SO_P1/10.0);
            P2 = float(QING_SO_P2/10.0);
        }
    }
}


void qing_mcost_to_disp::scanline_optimize_x(float *so_mcost_vol, float *in_mcost_vol, int xc, int r) {

    cout << "scanline optimize horizontal..\tr=" << r << endl;

    float * temp_mcost_vol = new float[m_total_size];
    if(0==temp_mcost_vol) {
        cerr << "bad alloc..." << endl;
        return ;
    }

    for(int d = 0; d < m_d; ++d) {
        for(int y = 0; y < m_h; ++y) {
            int idx = d * m_image_size + y * m_w + xc;
            temp_mcost_vol[idx] = in_mcost_vol[idx];
        }
    }

    for(int x = xc+r; 0<=x && m_w>x; x+=r) {
        for(int y = 0; y < m_h; ++y) {
            scanline_optim(temp_mcost_vol, in_mcost_vol, x, x-r, y, y);
        }
    }

    qing_add_so_mcost(so_mcost_vol, temp_mcost_vol, m_d, m_h, m_w);
    if(0!=temp_mcost_vol) {
        delete[] temp_mcost_vol; temp_mcost_vol = 0;
    }
}

void qing_mcost_to_disp::scanline_optimize_y(float * so_mcost_vol, float *in_mcost_vol, int yc, int r) {
    cout << "scanline optimize vertical..\tr=" << r << endl;

    float * temp_mcost_vol = new float[m_total_size];
    if(0==temp_mcost_vol) {
        cerr << "bad alloc..." << endl;
        return ;
    }

    for(int d = 0; d < m_d; ++d) {
        int idx_d = d * m_image_size + yc * m_w;
        float * ptr_yc_in = in_mcost_vol + idx_d;
        float * ptr_yc_temp = temp_mcost_vol + idx_d;
        memcpy(ptr_yc_temp, ptr_yc_in, sizeof(float)*m_w);
    }


    for(int y = yc+r; 0<=y && m_h>y; y+=r) {
        for(int x = 0; x < m_w; ++x) {
            scanline_optim(temp_mcost_vol, in_mcost_vol, x, x, y, y-r);
        }
    }

    qing_add_so_mcost(so_mcost_vol, temp_mcost_vol, m_d, m_h, m_w);
    if(0!=temp_mcost_vol) {
        delete[] temp_mcost_vol; temp_mcost_vol = 0;
    }
}

void qing_mcost_to_disp::scanline_optimize() {
    float * so_mcost = new float[m_total_size]; //scanline optimization matching cost
    std::fill(so_mcost, so_mcost+m_total_size, 0.f);

    scanline_optimize_x(so_mcost, m_filtered_mcost_l, 0, 1);                 //x-right
    scanline_optimize_x(so_mcost, m_filtered_mcost_l, m_w-1, -1);            //x-left
    scanline_optimize_y(so_mcost, m_filtered_mcost_l, 0, 1);                 //y-down
    scanline_optimize_y(so_mcost, m_filtered_mcost_l, m_h-1, -1);            //y-up

    std::fill(m_filtered_mcost_l, m_filtered_mcost_l + m_total_size, 0.f);
    memcpy(m_filtered_mcost_l, so_mcost, sizeof(float)*m_total_size);
}

