#include <string>
#include <vector>
#include <iostream>

#include "NPandSmoothingTools/SmoothingUtils.h"

int main (int argc, const char *argv[]) {
  std::vector<std::string> args;

  for (int i = 0; i < argc; ++i) args.emplace_back(argv[i]);

  if (args.size() < 2) {
    std::cerr << "This program needs the location of the input CDI file as \
it's first parameter. Aborting!" << std::endl;
    return 1;
  }

  auto        number_of_bins = 100;
  std::string container = ".*_SF",
              bNPset = "5,3,2",
              cNPset = "4,4,3",
              lightNPset = "12,5,4";
  float       b_bandwidth     = 0.4,
              c_bandwidth     = 0.4,
              light_bandwidth = 0.4;
  size_t      kernel_order   = 0;
  bool        sB             = true,
              sC             = true,
              sLight         = true;
  std::cout << "Running " << args[0] << " ..." << std::endl;

  if (args.size() >= 3) container = args[2];

  if (args.size() >= 4) number_of_bins = std::stoi(args[3]);

  if (args.size() >= 5) b_bandwidth = std::stof(args[4]);

  if (args.size() >= 6) c_bandwidth = std::stof(args[5]);

  if (args.size() >= 7) light_bandwidth = std::stof(args[6]);

  if (args.size() >= 8) kernel_order = std::stoi(args[7]);

  if (args.size() >= 9) sB = args[8] == "true";

  if (args.size() >= 10) bNPset = args[9];

  if (args.size() >= 11) sC = args[10] == "true";

  if (args.size() >= 12) cNPset = args[11];

  if (args.size() >= 13) sLight = args[12] == "true";

  if (args.size() >= 14) lightNPset = args[13];

  std::cout << "with the following parameters:" << std::endl;
  std::cout << "\tinput file - " << args[1] << std::endl;
  std::cout << "\tnumber of smoothed bins - " << number_of_bins << std::endl;
  std::cout << "\tkernel b-jet bandwidth - " << b_bandwidth << std::endl;
  std::cout << "\tkernel c-jet bandwidth - " << c_bandwidth << std::endl;
  std::cout << "\tkernel light-jet bandwidth - " << light_bandwidth << std::endl;
  std::cout << "\tkernel order - " << kernel_order << std::endl;

  if (sB) std::cout << "\twill smooth b-jet scale factors" << std::endl;

  if (sC) std::cout << "\twill smooth c-jet scale factors" << std::endl;

  if (sLight) std::cout << "\twill smooth light-jet scale factors" << std::endl;

  std::vector<int> bNPvec, cNPvec, lightNPvec;
  std::vector<std::string> tokens;
  tokens.push_back("");
  for (char c : bNPset) {
    if (c == ',') {
      tokens.push_back("");
      continue;
    }
    tokens.back().push_back(c);
  }
  for (auto &s : tokens) bNPvec.emplace_back(std::stoi(s));
  tokens.clear();
  tokens.push_back("");
  for (char c : cNPset) {
    if (c == ',') {
      tokens.push_back("");
      continue;
    }
    tokens.back().push_back(c);
  }
  for (auto &s : tokens) cNPvec.emplace_back(std::stoi(s));
  tokens.clear();
  tokens.push_back("");
  for (char c : lightNPset) {
    if (c == ',') {
      tokens.push_back("");
      continue;
    }
    tokens.back().push_back(c);
  }
  for (auto &s : tokens) lightNPvec.emplace_back(std::stoi(s));
  tokens.clear();

  Analysis::smoothContinuousCalibrations(args[1], container.c_str(), bNPvec, cNPvec, lightNPvec, number_of_bins, b_bandwidth, c_bandwidth, light_bandwidth, kernel_order, sB, sC, sLight);

  return 0;
} // main
