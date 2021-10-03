#include "profiledeserializer.hpp"
#include <fstream>
#include <string_view>
#include <string>

namespace ProfileDeserializer {

    std::vector<std::string_view> tokenize(const std::string_view& str, const char separator) {
        std::vector<std::string_view> tokens;
        size_t start;
        size_t end = 0;
   
        while ((start = str.find_first_not_of(separator, end)) != std::string::npos)
        {
            end = str.find(separator, start);
            tokens.push_back(str.substr(start, end - start));
        }
        return tokens;
    }

    std::optional<std::vector<Point>> open(const std::filesystem::path& path) {
        auto fp = std::ifstream(path.string());
        if (!fp.good()){
            printf("Could not open file. Error number: @TODO.");
            return {};
        }
        std::vector<Point> points;
        std::string string;
        //@TODO review
        while ((std::getline(fp, string)).good()) {
            std::string_view view(string);
            const auto tokens = tokenize(view, '\t');
            auto tokenIt = tokens.begin();
            //@TODO review
            while(tokenIt < tokens.cend()) {
                auto tokenSize = tokenIt->size();
                const double x = std::stod(tokenIt->data(), &tokenSize);
                ++tokenIt;
                tokenSize = tokenIt->size();
                const double y = std::stod(tokenIt->data(), &tokenSize);
                points.emplace_back(Point(x, y)); 
                ++tokenIt;
            }
        }
        return points;
    }
}
