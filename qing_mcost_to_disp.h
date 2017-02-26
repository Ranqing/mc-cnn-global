#ifndef QING_BIN_READER_H
#define QING_BIN_READER_H

//@ranqing
//@read matching cost files ".bin" to ".txt"

#include "../Qing/qing_common.h"
#include "../Qing/qing_string.h"
#include "../Qing/qing_file_reader.h"
#include "../Qing/qing_file_writer.h"

class qing_mcost_to_disp {
public :
    qing_mcost_to_disp(const int d, const int h, const int w) ;
    ~qing_mcost_to_disp() ;

    void read_from_mc_cnn_using_example_code(const string filename_l, const string filename_r);
    void read_from_mc_cnn(const string filename_l, const string filename_r );
    void mcost_to_disp();

private :
    float * m_mcost_l, * m_mcost_r;                            //matching cost
    int * m_disp_l, *m_disp_r, *m_disp;                               //disparity image
    int m_total_size, m_image_size;
    int m_d, m_h, m_w;                          //d: disparity_range, h: image height, w: image width

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
    m_disp_l = new int[m_image_size];
    m_disp_r = new int[m_image_size];
    if(0==m_disp_l || 0==m_disp_r) {
        cerr << "failed to initial disparity images." << endl;
        exit(-1);
    }
}

inline qing_mcost_to_disp::~qing_mcost_to_disp() {
    if(0!=m_mcost_l) { delete [] m_mcost_l;}
    if(0!=m_mcost_r) { delete [] m_mcost_r;}
    if(0!=m_disp_l)  { delete [] m_disp_l; }
    if(0!=m_disp_r)  { delete [] m_disp_r; }
}

#endif // QING_BIN_READER_H
