#include "json.h"

int main(int argc, char** argv) {
  vector<string> json_files;

  for (int i(1); i < argc; ++i) {
    json_files.push_back(argv[i]);
  }

  if (json_files.size() > 0) {
    cerr << "Previous output files specified." << endl;
    vector<DreamState> state(json_files.size(),20);

    for (size_t i(0); i < json_files.size(); ++i) {
      cerr << "Reading '" << json_files[i] << "'" << endl;
      string ll(last_line(json_files[i]));
      read_json(ll,state[i]);
      state[i].print();
    }

    for (size_t i(0); i < json_files.size(); ++i) {
      state[i].print();
    }
  }

  return 0;
}

