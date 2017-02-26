#include "qing_mcost_to_disp.h"

int main(int argc, char * argv[]) {

    const string filename_l = "../matching-cost/left.bin";
    const string filename_r = "../matching_cost/right.bin";
    int d = 70;
    int h = 370;
    int w = 1226;

    qing_mcost_to_disp solver(d, h, w);
    solver.read_from_mc_cnn_using_example_code(filename_l, filename_r);
    solver.mcost_to_disp();
    return 1;
}
