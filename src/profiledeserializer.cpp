#include "profiledeserializer.hpp"
#include <string_view>
#include <stdlib.h>
#include <errno.h>

namespace ProfileDeserializer {

    std::vector<std::string_view> tokenize(const std::string_view& inString, const char separator) {
        std::vector<std::string_view> tokens;
        auto string = inString;
        auto pos = string.find(separator);
        while(pos != std::string::npos) {
            tokens.emplace_back(string.substr(0, pos));
            string.remove_prefix(pos + 1);
            pos = string.find(separator);
            if (pos = std::string::npos) {
                pos = string.find('\n');
            }
        }
        return tokens;
    }

    std::optional<std::vector<Point>> open(const std::filesystem::path& path) {
        errno = 0;
        FILE* fp = fopen(path.c_str(), "r");
        if (fp == NULL){
            auto er = errno;
            printf("Could not open file. Error number: %i.", er);
            return {};
        }
        char* line = nullptr;
        size_t len = 0;
        char * tokens = nullptr;
        std::vector<Point> points;
        while ((getline(&line, &len, fp)) != -1) {
            std::string_view view(line);
            const auto tokens = tokenize(view, '\t');
            auto tokenIt = tokens.begin();
            while(tokenIt < tokens.cend()) {
                auto tokenSize = tokenIt->size();
                double x = std::stod(tokenIt->data(), &tokenSize);
                ++tokenIt;
                tokenSize = tokenIt->size();
                double y = std::stod(tokenIt->data(), &tokenSize);
                points.emplace_back(Point(x, y)); 
                ++tokenIt;
            }
        }
        if (line) {
            free(line);
        }
        fclose(fp);
        return points;
    }
}
