/// @brief Quick smoke-test for filter_design_construct.

#include <iostream>
#include <multiplierless/lowpass_oracle.hpp>

/**
 * @brief Verify that filter_design_construct produces valid constraint matrices.
 * @return 0 on success.
 */
int main() {
    std::cout << "Creating Fdc..." << std::endl;
    auto Fdc = filter_design_construct(32);
    std::cout << "Fdc.Ap.rows=" << Fdc.Ap.rows() << " cols=" << Fdc.Ap.cols() << std::endl;
    std::cout << "Fdc.As.rows=" << Fdc.As.rows() << " cols=" << Fdc.As.cols() << std::endl;
    std::cout << "Fdc.Anr.rows=" << Fdc.Anr.rows() << " cols=" << Fdc.Anr.cols() << std::endl;
    std::cout << "OK" << std::endl;
    return 0;
}
