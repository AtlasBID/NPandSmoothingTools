#include <string>
#include <vector>
#include <iostream>

#include "NPandSmoothingTools/SmoothingUtils.h"

int main (int argc, const char *argv[]) {
  std::vector<std::string> args;
  for (auto i = 0; i < argc; ++i)
    args.emplace_back(argv[i]);

  auto merge_num = -1;
  auto preserve_error = false;
  std::cout << "Running " << args[0] << " ..." << std::endl;

  if (args.size() < 4) {
    std::cerr << "This program needs the location of the input CDI file as "
      "it's first parameter, an unsmoothed (reference) CDI file for the "
      "second parameter, and the CDI object to analyze as the third. "
      "Aborting!" << std::endl;
    return 1;
  }
  if (args.size() > 4)
    merge_num = std::stoi(args[4]);
  if (args.size() > 5)
    preserve_error = (args[5] == "true");

  // std::cout << "with the following parameters:" << std::endl;
  // std::cout << "\tinput file - " << args[1] << std::endl;
  // std::cout << "\tnumber of smoothed bins - " << number_of_bins << std::endl;
  // std::cout << "\tkernel bandwidth - " << bandwidth << std::endl;

  Analysis::eigenDecompCalib(args[1], args[2], args[3], merge_num, preserve_error);
  //\"MV1/AntiKt4TopoLCJVF0_5/0_7892/B/combined_pdf_dijet_10_SF\")");
  //\"MV1/AntiKt4TopoLCJVF0_5/0_3511/B/KinSel_dilep_SF\")");
  //\"MV1/AntiKt4TopoLCJVF0_5/0_3511/C/DStar_ttbar_PDF_7b_rebin_SF\")");

  return 0;
}
