#pragma once
#include <vector>
#include <filesystem>
#include <optional>
#include "point.hpp"

namespace ProfileDeserializer {
    std::optional<std::vector<Point>> open(const std::filesystem::path& path);
}
