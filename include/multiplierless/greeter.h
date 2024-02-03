#pragma once

#include <string>

namespace multiplierless {

    /**  Language codes to be used with the Multiplierless class */
    enum class LanguageCode { EN, DE, ES, FR };

    /**
     * @brief A class for saying hello in multiple languages
     */
    class Multiplierless {
        std::string name;

      public:
        /**
         * @brief Creates a new multiplierless
         * @param[in] name the name to greet
         */
        Multiplierless(std::string name);

        /**
         * @brief Creates a localized string containing the greeting
         * @param[in] lang the language to greet in
         * @return a string containing the greeting
         */
        std::string greet(LanguageCode lang = LanguageCode::EN) const;
    };

}  // namespace multiplierless
