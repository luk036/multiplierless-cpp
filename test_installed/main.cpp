#include <multiplierless/version.h>

#include <iostream>

auto main() -> int {
    const auto ok = (MULTIPLIERLESS_VERSION_MAJOR >= 1);
    std::cout << "multiplierless installed test: version " << MULTIPLIERLESS_VERSION << "\n";
    return ok ? 0 : 1;
}
