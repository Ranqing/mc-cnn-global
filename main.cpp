#include "qing_mcost_to_disp.h"

int main(int argc, char * argv[]) {

//    const string filename_l = "../matching-cost-half/left.bin";
//    const string filename_r = "../matching-cost-half/right.bin";
//    const string imagename_l = "../input/half_kittiL.png";
//    const string imagename_r = "../input/half_kittiR.png";

//    int d = 40;
//    int h = 185;
//    int w = 613;
//    float sigma_range = 0.06;
//    float sigma_spatial = 0.03;
//    int wnd = 21;

    const string filename_l = "../matching-cost/left.bin";
    const string filename_r = "../matching-cost/right.bin";
    const string imagename_l = "../input/kittiL.png";
    const string imagename_r = "../input/kittiR.png";

    int d = 70;
    int h = 370;
    int w = 1226;
    float sigma_range = 0.06;
    float sigma_spatial = 0.03;
    int wnd = 21;

    qing_mcost_to_disp solver(d, h, w);
    QingTimer timer;
    solver.read_image(imagename_l, imagename_r);
    cout << "loading image.." << timer.duration()*1000 << " ms" << endl;

# if 0
    timer.restart();
    solver.read_from_mc_cnn_using_example_code(filename_l, filename_r);
    solver.remove_mcost_nan();
    cout << "loading matching cost.." << timer.duration()*1000 << " ms" << endl;

    solver.get_weighted_table(sigma_range, sigma_spatial);
    cout << "weighted table calculation done..." << endl;
# endif

# if 0
    timer.restart();
    solver.mcost_aggregation(wnd);
    cout << "matching cost aggregation.." << timer.duration()*1000 << " ms" << endl;
# endif

# if 0
    timer.restart();
    solver.directional_mcost_aggregation(wnd);
    cout << "directional matching cost aggregation.." << timer.duration()*1000 << " ms" << endl;
    qing_create_dir("qing-matching-cost-half");
    solver.save_filtered_mcost("qing-matching-cost-half");
# endif

# if 1
    string mcost_folder = "../qing-matching-cost";  //../qing-matching-cost-half
    solver.read_filtered_mcost(mcost_folder);
    solver.scanline_optimize();
   // solver.semi_global();
# endif

    timer.restart();
    solver.mcost_to_disp(255/d, "d_bf_disp_so.png");
    cout << "matching cost to disparity.." << timer.duration()*1000 << " ms" << endl;
    return 1;
}
