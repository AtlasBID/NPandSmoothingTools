
#ifndef TRAITS_HELPERS_H
#define TRAITS_HELPERS_H

#include <type_traits>

namespace Analysis {

// variatic template type_traits helpers

template <typename toComp, typename... Args>
struct are_convertible;

template <typename toComp, typename Arg>
struct are_convertible <toComp, Arg>
{
  static const auto value = std::is_convertible<toComp, Arg>::value;
};
template <typename toComp, typename Arg, typename... Args>
struct are_convertible <toComp, Arg, Args...>
{
  static const auto value = std::is_convertible<toComp, Arg>::value && are_convertible<toComp, Args...>::value;
};


template <typename toComp, typename... Args>
struct are_assignable;

template <typename toComp, typename Arg>
struct are_assignable <toComp, Arg>
{
  static const auto value = std::is_assignable<toComp, Arg>::value;
};
template <typename toComp, typename Arg, typename... Args>
struct are_assignable <toComp, Arg, Args...>
{
  static const auto value = std::is_assignable<toComp, Arg>::value && are_assignable<toComp, Args...>::value;
};


template <typename... Args>
struct are_arithmetic;

template <typename Arg>
struct are_arithmetic <Arg>
{
  static const auto value = std::is_arithmetic<Arg>::value;
};
template <typename Arg, typename... Args>
struct are_arithmetic<Arg, Args...>
{
  static const auto value = std::is_arithmetic<Arg>::value && are_arithmetic<Args...>::value;
};

template <typename... Args>
struct are_integral;

template <typename Arg>
struct are_integral <Arg>
{
  static const auto value = std::is_integral<Arg>::value;
};
template <typename Arg, typename... Args>
struct are_integral<Arg, Args...>
{
  static const auto value = std::is_integral<Arg>::value && are_integral<Args...>::value;
};

}

#endif // include gaurds
