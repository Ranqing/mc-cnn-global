#include "qing_mcost_to_disp.h"

void qing_mcost_to_disp::read_from_mc_cnn(const string filename_l, const string filename_r) {
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

# if 1
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

void qing_mcost_to_disp::read_from_mc_cnn_using_example_code(const string filename_l, const string filename_r) {
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

# if 1
     string out_file_l = getFilePrefix(filename_l) + "_example.txt";
     string out_file_r = getFilePrefix(filename_r) + "_example.txt";
     qing_write_txt(out_file_l, m_mcost_l, m_total_size, m_h * m_w);
     qing_write_txt(out_file_r, m_mcost_r, m_total_size, m_h * m_w);
# endif
}

void qing_mcost_to_disp::mcost_to_disp() {
}

