#ifndef OPT_PARSE_HPP_
#define OPT_PARSE_HPP_

#include <cstdio>
#include <unordered_map>
#include <stdexcept>
using namespace std;

class OptionParser {

private:
    unordered_map<string, string> optionsMap;

public:
    OptionParser(int argc, char** argv) {
        if (argc % 2 == 0) {
            throw invalid_argument("Expected a set of options of the form (key, value)");
        }

        for (int i = 1; i < argc; i += 2) {
            optionsMap[argv[i]] = argv[i + 1];
        }
    }

    inline void ensureOptionExists(const string& option) {
        if (optionsMap.find(option) == optionsMap.end()) {
            string message = "Missing the " + option + " parameter";
            throw invalid_argument(message);
            
        }
    }

    string& getOptionValue(const string& option) {
        ensureOptionExists(option);
        return optionsMap[option];
    }
};

#endif