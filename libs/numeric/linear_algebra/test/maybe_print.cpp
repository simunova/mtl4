#include <iostream>

// From Doug

auto concept Any<typename T> { }

auto concept Streamable<typename T> {
  std::ostream& operator<<(std::ostream &, const T&);
}

template<Any T> void print(const T&) {
  std::cout << "<not printable>" << std::endl;
}

template<Any T> requires Streamable<T> void print(const T& x) {
  std::cout << x << std::endl;
}

template<Any T>
void maybe_print(const T& x) {
  print(x);
}

struct NotPrintable { };

int main() {
  print(17);
  print(NotPrintable());
  return 0;
}
