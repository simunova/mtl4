template <typename Op, typename T, T b>
struct op_right_constant
{
	static apply(T a) { return Op()(a, b); }
};

template <typename Op, typename T>
struct op_right_constant<T, a, static_identity<Op, T>::value>
{	
	static apply(T a) {return a; }
}

struct add
{
	template <typename T>
	static apply(T a, T b) { return a + b; }
};

template <typename Op, typename T>
struct static_identity {};

template <typename T>
struct static_identity<add, T>
{
	const static T value= 0;
};

template <>
struct static_identity<add, std::string>
{
	const static std::string value;
};

static_identity<add, std::string>::value("");

template <typename T, T y>
inline add_f(T x)
{
		return op_right_constant<add, y>::apply(x);
}

template <typename T>
int f(T);

template <typename T> requires C1<T>
int f(T);

template <typename T> requires C2<T>
int f(T);

template <typename T> requires C3<T>
int f(T);

template <typename T> requires C4<T>
int f(T);



template <typename T>
typename enable_if<mpl::bool<!C1<T>::value && C2<T>::value>, int>::type
f(T);









int main()
{
	string x("text");
	const string y("");
	string z= op_right_constant<add, y>::apply(x);
	string a= add_f<string, y>(x);
}
		