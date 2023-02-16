#include <multiplierless/greeter.h> // for Multiplierless, multiplierless

#include <iosfwd>      // for string
#include <string>      // for basic_string
#include <type_traits> // for move

using namespace multiplierless;

Multiplierless::Multiplierless(std::string _name) : name(std::move(_name)) {}

// std::string Multiplierless::greet(LanguageCode lang) const {
//     switch (lang) {
//         default:
//         case LanguageCode::EN:
//             return fmt::format("Hello, {}!", name);
//         case LanguageCode::DE:
//             return fmt::format("Hallo {}!", name);
//         case LanguageCode::ES:
//             return fmt::format("¡Hola {}!", name);
//         case LanguageCode::FR:
//             return fmt::format("Bonjour {}!", name);
//     }
// }
