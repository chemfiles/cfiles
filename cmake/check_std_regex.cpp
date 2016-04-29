#include <regex>

int main() {
    auto words = std::string("Quick brown fox.");

    std::regex words_regex("[^\\s]+");
    auto words_begin = std::sregex_iterator(words.begin(), words.end(), words_regex);
    auto words_end = std::sregex_iterator();

    return std::distance(words_begin, words_end);
}
