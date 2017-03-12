#include "qing_mcost_to_disp.h"

int main(int argc, char * argv[]) {

    //    const string filename_l = "../matching-cost-half/left.bin";
    //    const string filename_r = "../matching-cost-half/right.bin";
    //    const string imagename_l = "../input/half_kittiL.png";
    //    const string imagename_r = "../input/half_kittiR.png";

    //    int d = 40;
    //    int h = 185;
    //    int w = 613;
    //    float sigma_range = 0.08;
    //    float sigma_spatial = 0.03;

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

    timer.restart();
    solver.read_from_mc_cnn_using_example_code(filename_l, filename_r);
    // solver.read_from_mc_cnn(filename_l, filename_r);

    solver.remove_mcost_nan();
    cout << "loading matching cost.." << timer.duration()*1000 << " ms" << endl;

# if 0
    timer.restart();
    solver.mcost_aggregation(sigma_range, sigma_spatial);
    cout << "matching cost aggregation.." << timer.duration()*1000 << " ms" << endl;
# endif

# if 1
    timer.restart();
    solver.directional_mcost_aggregation(sigma_range, sigma_spatial, wnd);
    cout << "directional matching cost aggregation.." << timer.duration()*1000 << " ms" << endl;
# endif
# if 1
    timer.restart();
    const int scale = 255/d;
    solver.mcost_to_disp(scale, "_disp.png");
    cout << "matching cost to disparity.." << timer.duration()*1000 << " ms" << endl;
# endif
    return 1;
}
