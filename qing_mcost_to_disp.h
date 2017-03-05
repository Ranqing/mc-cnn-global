#ifndef QING_BIN_READER_H
#define QING_BIN_READER_H

//@ranqing

#include "../../../Qing/qing_common.h"
#include "../../../Qing/qing_string.h"
#include "../../../Qing/qing_timer.h"
#include "../../../Qing/qing_dir.h"
#include "../../../Qing/qing_file_reader.h"
#include "../../../Qing/qing_file_writer.h"
#include "../../../Qing/qing_check.h"
#include "../../../Qing/qing_macros.h"
#include "../../../Qing/qing_matching_cost.h"
#include "../../../Qing/qing_bilateral_filter.h"

class qing_mcost_to_disp {
public :
    qing_mcost_to_disp(const int d, const int h, const int w) ;
    ~qing_mcost_to_disp() ;

    void read_image(const string filename_l, const string filename_r);
    void read_from_mc_cnn_using_example_code(const string filename_l, const string filename_r);
    void read_from_mc_cnn(const string filename_l, const string filename_r );
    void mcost_to_disp(const int scale, const string savename);

    void mcost_aggregation(float sigma_range, float sigma_spatial) ;
    void adaptive_weight_filter(float *out, float *in, uchar * gray_l, uchar * gray_r, int d, int wnd);
    //void adaptive_weight_filter(float * out, float * in, uchar * colors_l, uchar * colors_r,  int d, int h, int w, int wnd, float sigma_spatial, float sigma_range );


    void directional_mcost_aggregation(float sigma_range, float sigma_spatial);
    void directional_adaptive_weight_filter(float *out, float *in, uchar * gray_l, uchar * gray_r, int d, int wnd);

private:
    float * m_mcost_l, * m_mcost_r;                                    //matching cost
    float * m_filtered_mcost_l, * m_filtered_mcost_r;                  //filtered matching cost
    uchar * m_image_l, * m_image_r;                                    //left and right image
    int * m_disp_l, * m_disp_r, * m_disp;                              //disparity image
    int m_total_size, m_image_size;                                    //total_size = d * h * w, image_size = h * w
    int m_d, m_h, m_w, m_channels;                                     //d: disparity_range, h: image height, w: image width

    float * m_range_table;
    float * m_spatial_table;

    Mat m_mat_l, m_mat_r;
};

inline qing_mcost_to_disp::qing_mcost_to_disp(const int d, const int h, const int w): m_d(d), m_h(h), m_w(w) {
    m_total_size = m_d * m_h * m_w;
    m_image_size = m_h * m_w;
    m_mcost_l = new float[m_total_size];
    m_mcost_r = new float[m_total_size];
    if(0==m_mcost_l || 0==m_mcost_r) {
        cerr << "failed to initial matching cost array." << endl;
        exit(-1);
    }
//initialization in aggregation
//    m_filtered_mcost_l = new float[m_total_size];
//    m_filtered_mcost_r = new float[m_total_size];
//    if(0==m_filtered_mcost_l || 0==m_filtered_mcost_r) {
//        cerr << "failed to initial filtered matching cost array." << endl;
//        exit(-1);
//    }
    m_filtered_mcost_l = 0;
    m_filtered_mcost_r = 0;

    m_disp_l = new int[m_image_size];
    m_disp_r = new int[m_image_size];
    m_disp = new int[m_image_size];
    if(0==m_disp_l || 0==m_disp_r) {
        cerr << "failed to initial disparity images." << endl;
        exit(-1);
    }
}

inline qing_mcost_to_disp::~qing_mcost_to_disp() {
//    if(0!=m_mcost_l) { delete [] m_mcost_l;}
//    if(0!=m_mcost_r) { delete [] m_mcost_r;}
//    if(0!=m_disp_l)  { delete [] m_disp_l; }
//    if(0!=m_disp_r)  { delete [] m_disp_r; }
}


inline void qing_mcost_to_disp::read_image(const string filename_l, const string filename_r) {
    m_mat_l = imread(filename_l, CV_LOAD_IMAGE_UNCHANGED);
    if(0==m_mat_l.data) {
        cerr << "failed to open " << filename_l << endl;
        return ;
    }
    m_mat_r = imread(filename_r, CV_LOAD_IMAGE_UNCHANGED);
    if(0==m_mat_r.data) {
        cerr << "failed to open " << filename_r << endl;
        return ;
    }
    m_channels = m_mat_l.channels();
    cout << m_mat_l.channels() << endl;

    m_image_l = new uchar[m_image_size * m_channels];
    m_image_r = new uchar[m_image_size * m_channels];
    if(0==m_image_l || 0==m_image_r) {
        cerr << "failed to initial bgr images." << endl;
        exit(-1);
    }
    memcpy(m_image_l, m_mat_l.data, sizeof(uchar) * m_image_size * m_channels);
    memcpy(m_image_r, m_mat_r.data, sizeof(uchar) * m_image_size * m_channels);

# if 1
   //cout << "COOL" << endl;
   Mat test_mat_l(m_h, m_w, m_mat_l.type());
   Mat test_mat_r(m_h, m_w, m_mat_r.type());
   memcpy(test_mat_l.data, m_image_l, sizeof(uchar) * m_image_size * m_channels);
   memcpy(test_mat_r.data, m_image_r, sizeof(uchar) * m_image_size * m_channels);
   imshow("test_mat_l", test_mat_l);
   imshow("test_mat_r", test_mat_r);
   waitKey(0);
   destroyAllWindows();
# endif
}


//@read matching cost files ".bin" to ".txt"
inline void qing_mcost_to_disp::read_from_mc_cnn(const string filename_l, const string filename_r) {
    if( false == qing_check_file_suffix(filename_l, "bin") ) {
        cerr << "failed to open " << filename_l << endl;
        return ;
    }
    if( false == qing_check_file_suffix(filename_r, "bin") ) {
        cerr << "failed to open " << filename_r << endl;
        return ;
    }

    qing_read_bin(filename_l, m_mcost_l, m_total_size) ;
    qing_read_bin(filename_r, m_mcost_r, m_total_size) ;

# if 0
    ////    string prefix = qing_get_file_prefix(filename);
    ////    string out_filename;
    ////    for(int i = 0; i < m_d; ++i) {
    ////        out_filename = prefix + "_mc_" + qing_int_2_format_string(i,2,'0') + ".txt";
    ////        qing_write_txt(out_filename, m_mcost + i * m_image_size, m_image_size, m_w);
    ////    }

    string out_file_l = getFilePrefix(filename_l) + ".txt";
    string out_file_r = getFilePrefix(filename_r) + ".txt";
    qing_write_txt(out_file_l, m_mcost_l, m_total_size, m_h * m_w);
    qing_write_txt(out_file_r, m_mcost_r, m_total_size, m_h * m_w);
# endif
}

inline void qing_mcost_to_disp::read_from_mc_cnn_using_example_code(const string filename_l, const string filename_r) {
    int fd_l = open(filename_l.c_str(), O_RDONLY);
    if(-1==fd_l) {
        cerr << "failed to open " << filename_l;
        return ;
    }
    m_mcost_l = (float *)mmap(NULL, 1*m_d*m_w*m_h*sizeof(float), PROT_READ, MAP_SHARED, fd_l, 0);    //fd: file descriptor
    close(fd_l);

    int fd_r = open(filename_r.c_str(), O_RDONLY);
    if(-1==fd_r) {
        cerr << "failed to open " << filename_r;
        return ;
    }
    m_mcost_r = (float *)mmap(NULL, 1*m_d*m_w*m_h*sizeof(float), PROT_READ, MAP_SHARED, fd_r, 0);
    close(fd_r);

# if 0
    qing_create_dir("matching_cost");
    for(int d = 0; d < m_d; ++d) {
        float * mcost = m_mcost_l + d * m_image_size;
        string out_file = "matching_cost/mcost_l_" + qing_int_2_string(d) + ".jpg";
        qing_save_mcost_jpg(out_file, mcost, m_w, m_h);
    }
# endif
}

#endif // QING_BIN_READER_H
