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

    const string mcost_folder = "../matching-cost-cnn/";
    const string image_folder = "../input/";


//    int d = 70;
//    int h = 370;
//    int w = 1226;
    int d = 255;
    int h = 375;
    int w = 1242;
    float sigma_range = 0.06;
    float sigma_spatial = 0.03;
    int wnd = 21;

    qing_mcost_to_disp solver(d, h, w);

    QingTimer timer;
    string imagename_l = "kittiL_2.png";
    string imagename_r = "kittiR_2.png";
    solver.read_image(image_folder + imagename_l, image_folder + imagename_r);
    cout << "loading image.." << timer.duration()*1000 << " ms" << endl;

//# if 1
//    timer.restart();
//    string mcostname_l = "left.bin";
//    string mcostname_r = "right.bin";
//    solver.read_from_mc_cnn_using_example_code(mcost_folder + mcostname_l, mcost_folder + mcostname_r);
//    solver.remove_mcost_nan();

//    //shift image up in y-direction
//    mcostname_l = "left_u.bin";
//    mcostname_r = "right_u.bin";
//    cout << mcostname_l << '\t' << mcostname_r << endl;
//    solver.read_from_mc_cnn_cmp(mcost_folder + mcostname_l, mcost_folder + mcostname_r);     //compare each loading matching cost

//    //shift image down in y-direction
//    mcostname_l = "left_d.bin";
//    mcostname_r = "right_d.bin";
//    cout << mcostname_l << '\t' << mcostname_r << endl;
//    solver.read_from_mc_cnn_cmp(mcost_folder + mcostname_l, mcost_folder + mcostname_r);    //compare each loading matching cost
//    cout << "loading matching cost.." << timer.duration()*1000 << " ms" << endl;

//    solver.get_weighted_table(sigma_range, sigma_spatial);
//    cout << "weighted table calculation done..." << endl;
//# endif

# if 1
    timer.restart();
    string mcostname_l = "left_2.bin";
    string mcostname_r = "right_2.bin";
    solver.read_from_mc_cnn_using_example_code(mcost_folder + mcostname_l, mcost_folder + mcostname_r);
    solver.remove_mcost_nan();
    cout << "loading matching cost.." << timer.duration()*1000 << " ms" << endl;

    solver.get_weighted_table(sigma_range, sigma_spatial);
    cout << "weighted table calculation done..." << endl;
# endif

# if 1
    timer.restart();
    solver.mcost_aggregation(wnd);
    cout << "matching cost aggregation.." << timer.duration()*1000 << " ms" << endl;
# endif

# if 1
    timer.restart();
    solver.directional_mcost_aggregation(wnd);
    cout << "directional matching cost aggregation.." << timer.duration()*1000 << " ms" << endl;
    qing_create_dir("qing-matching-cost");
    solver.save_filtered_mcost("qing-matching-cost-2");
# endif

# if 1
   // string mcost_folder = "../qing-matching-cost";  //../qing-matching-cost-half
   // solver.read_filtered_mcost(mcost_folder);
    solver.scanline_optimize();
    qing_create_dir("qing-matching-cost-so-rect");
    solver.save_filtered_mcost("qing-matching-cost-so-rect");
   // solver.semi_global();
# endif

    timer.restart();
    solver.mcost_to_disp(255/d, "disp.png");
    cout << "matching cost to disparity.." << timer.duration()*1000 << " ms" << endl;

    //solver.evaluate();
    return 1;
}
