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
  float       bandwidth      = 0.4;
  size_t      kernel_order   = 0;
  bool        sB             = true,
              sC             = true,
              sLight         = true;
  std::cout << "Running " << args[0] << " ..." << std::endl;

  if (args.size() >= 3) container = args[2];

  if (args.size() >= 4) number_of_bins = std::stoi(args[3]);

  if (args.size() >= 5) bandwidth = std::stof(args[4]);

  if (args.size() >= 6) kernel_order = std::stoi(args[5]);

  if (args.size() >= 7) sB = args[6] == "true";

  if (args.size() >= 8) bNPset = args[7];

  if (args.size() >= 9) sC = args[8] == "true";

  if (args.size() >= 10) cNPset = args[9];

  if (args.size() >= 11) sLight = args[10] == "true";

  if (args.size() >= 12) lightNPset = args[11];

  std::cout << "with the following parameters:" << std::endl;
  std::cout << "\tinput file - " << args[1] << std::endl;
  std::cout << "\tnumber of smoothed bins - " << number_of_bins << std::endl;
  std::cout << "\tkernel bandwidth - " << bandwidth << std::endl;
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

  Analysis::smoothContinuousCalibrations(args[1], container.c_str(), bNPvec, cNPvec, lightNPvec, number_of_bins, bandwidth, kernel_order, sB, sC, sLight);

  return 0;
} // main
