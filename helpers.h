#pragma once

#ifndef MINIZ_NO_ZLIB_COMPATIBLE_NAMES
#define MINIZ_NO_ZLIB_COMPATIBLE_NAMES
#endif
#include "basisu_miniz.h"

#include "lodepng.h"
#include "lodepng_util.h"

#ifndef TINYEXR_USE_ZFP
#define TINYEXR_USE_ZFP (1)
#endif
#include "tinyexr.h"

#ifdef _MSC_VER
#pragma warning (disable: 4201) // nonstandard extension used: nameless struct/union
#pragma warning (disable: 4127) // conditional expression is constant
#endif

namespace basisu 
{
	static inline uint16_t byteswap16(uint16_t x) { return static_cast<uint16_t>((x << 8) | (x >> 8)); }
	static inline uint32_t byteswap32(uint32_t x) { return ((x << 24) | ((x << 8) & 0x00FF0000) | ((x >> 8) & 0x0000FF00) | (x >> 24)); }

	typedef std::vector<uint8_t> uint8_vec;
	typedef std::vector<float> float_vec;
	typedef std::vector<uint16_t> uint16_vec;
	typedef std::vector<uint32_t> uint32_vec;

	inline int posmod(int x, int y)
	{
		if (x >= 0)
			return (x < y) ? x : (x % y);
		int m = (-x) % y;
		return (m != 0) ? (y - m) : m;
	}

	inline uint32_t clz(uint32_t x)
	{
		if (!x)
			return 32;

		uint32_t n = 0;
		while ((x & 0x80000000) == 0)
		{
			x <<= 1u;
			n++;
		}

		return n;
	}

	inline uint32_t rot_left(uint32_t v, uint32_t s)
	{
		assert(s <= 32);
		if ((s) && (s != 32))
			v = (v << s) | (v >> (32 - s));
		return v;
	}

	inline uint32_t rot_right(uint32_t v, uint32_t s)
	{
		assert(s <= 32);
		if ((s) && (s != 32))
			v = (v >> s) | (v << (32 - s));
		return v;
	}

	inline uint16_t rot_left16(uint16_t v, uint16_t s)
	{
		assert(s <= 16);
		if ((s) && (s != 16))
			v = (v << s) | (v >> (16 - s));
		return v;
	}

	inline uint16_t rot_right16(uint16_t v, uint16_t s)
	{
		assert(s <= 16);
		if ((s) && (s != 16))
			v = (v >> s) | (v << (16 - s));
		return v;
	}

	inline uint8_t clamp255(int32_t i)
	{
		return (uint8_t)((i & 0xFFFFFF00U) ? (~(i >> 31)) : i);
	}

	enum eZero
	{
		cZero = 0
	};

	enum eNoClamp
	{
		cNoClamp = 0
	};

	template <typename T> inline void clear_obj(T& obj) { memset(&obj, 0, sizeof(obj)); }

	template <typename T0, typename T1> inline T0 lerp(T0 a, T0 b, T1 c) { return a + (b - a) * c; }

	template <typename S> inline S maximum(S a, S b) { return (a > b) ? a : b; }
	template <typename S> inline S maximum(S a, S b, S c) { return maximum(maximum(a, b), c); }
	template <typename S> inline S maximum(S a, S b, S c, S d) { return maximum(maximum(maximum(a, b), c), d); }

	template <typename S> inline S minimum(S a, S b) { return (a < b) ? a : b; }
	template <typename S> inline S minimum(S a, S b, S c) { return minimum(minimum(a, b), c); }
	template <typename S> inline S minimum(S a, S b, S c, S d) { return minimum(minimum(minimum(a, b), c), d); }

	inline float clampf(float value, float low, float high) { if (value < low) value = low; else if (value > high) value = high;	return value; }
	inline float saturate(float value) { return clampf(value, 0, 1.0f); }
	inline uint8_t minimumub(uint8_t a, uint8_t b) { return (a < b) ? a : b; }
	inline uint32_t minimumu(uint32_t a, uint32_t b) { return (a < b) ? a : b; }
	inline int32_t minimumi(int32_t a, int32_t b) { return (a < b) ? a : b; }
	inline float minimumf(float a, float b) { return (a < b) ? a : b; }
	inline uint8_t maximumub(uint8_t a, uint8_t b) { return (a > b) ? a : b; }
	inline uint32_t maximumu(uint32_t a, uint32_t b) { return (a > b) ? a : b; }
	inline int32_t maximumi(int32_t a, int32_t b) { return (a > b) ? a : b; }
	inline float maximumf(float a, float b) { return (a > b) ? a : b; }
	inline int squarei(int i) { return i * i; }
	inline float squaref(float i) { return i * i; }
	template<typename T> inline T square(T a) { return a * a; }

	template <typename S> inline S clamp(S value, S low, S high) { return (value < low) ? low : ((value > high) ? high : value); }

	template<typename T> inline void clear_vector(T& vec) { vec.erase(vec.begin(), vec.end()); }
	template<typename T> inline typename T::value_type* enlarge_vector(T& vec, size_t n) { size_t cs = vec.size(); vec.resize(cs + n); return &vec[cs]; }

	typedef uint16_t half_float;

	const double MIN_DENORM_HALF_FLOAT = 0.000000059604645; // smallest positive subnormal number
	const double MIN_HALF_FLOAT = 0.00006103515625; // smallest positive normal number
	const double MAX_HALF_FLOAT = 65504.0; // largest normal number

	inline bool in_range(int v, int l, int h)
	{
		return (v >= l) && (v <= h);
	}

	inline uint32_t get_bits(uint32_t val, int low, int high)
	{
		const int num_bits = (high - low) + 1;
		assert(in_range(num_bits, 1, 32));

		val >>= low;
		if (num_bits != 32)
			val &= ((1u << num_bits) - 1);

		return val;
	}

	inline bool is_half_inf_or_nan(half_float v)
	{
		return get_bits(v, 10, 14) == 31;
	}

	inline bool is_half_denorm(half_float v)
	{
		int e = (v >> 10) & 31;
		return !e;
	}

	inline int get_half_exp(half_float v)
	{
		int e = ((v >> 10) & 31);
		return e ? (e - 15) : -14;
	}

	inline int get_half_mantissa(half_float v)
	{
		if (is_half_denorm(v))
			return v & 0x3FF;
		return (v & 0x3FF) | 0x400;
	}

	inline float get_half_mantissaf(half_float v)
	{
		return ((float)get_half_mantissa(v)) / 1024.0f;
	}

	inline int get_half_sign(half_float v)
	{
		return v ? ((v & 0x8000) ? -1 : 1) : 0;
	}

	half_float float_to_half(float val)
	{
		union { float f; int32_t i; uint32_t u; } fi = { val };
		const int flt_m = fi.i & 0x7FFFFF, flt_e = (fi.i >> 23) & 0xFF, flt_s = (fi.i >> 31) & 0x1;
		int s = flt_s, e = 0, m = 0;

		// inf/NaN
		if (flt_e == 0xff)
		{
			e = 31;
			if (flt_m != 0) // NaN
				m = 1;
		}
		// not zero or denormal
		else if (flt_e != 0)
		{
			int new_exp = flt_e - 127;
			if (new_exp > 15)
				e = 31;
			else if (new_exp < -14)
				m = lrintf((1 << 24) * fabsf(fi.f));
			else
			{
				e = new_exp + 15;
				m = lrintf(flt_m / (float)(1 << 13));
			}
		}

		assert(0 <= m && m <= 1024);
		if (m == 1024)
		{
			e++;
			m = 0;
		}

		assert((s >= 0) && (s <= 1));
		assert((e >= 0) && (e <= 31));
		assert((m >= 0) && (m <= 1023));

		half_float result = (half_float)((s << 15) | (e << 10) | m);
		return result;
	}

	inline float half_to_float(half_float hval)
	{
		union { float f; uint32_t u; } x = { 0 };

		uint32_t s = ((uint32_t)hval >> 15) & 1;
		uint32_t e = ((uint32_t)hval >> 10) & 0x1F;
		uint32_t m = (uint32_t)hval & 0x3FF;

		if (!e)
		{
			if (!m)
			{
				// +- 0
				x.u = s << 31;
				return x.f;
			}
			else
			{
				// denormalized
				while (!(m & 0x00000400))
				{
					m <<= 1;
					--e;
				}

				++e;
				m &= ~0x00000400;
			}
		}
		else if (e == 31)
		{
			if (m == 0)
			{
				// +/- INF
				x.u = (s << 31) | 0x7f800000;
				return x.f;
			}
			else
			{
				// +/- NaN
				x.u = (s << 31) | 0x7f800000 | (m << 13);
				return x.f;
			}
		}

		e = e + (127 - 15);
		m = m << 13;

		assert(s <= 1);
		assert(m <= 0x7FFFFF);
		assert(e <= 255);

		x.u = m | (e << 23) | (s << 31);
		return x.f;
	}

	template <uint32_t N, typename T>
	class vec
	{
	protected:
		T m_v[N];

	public:
		enum { num_elements = N };
		typedef T scalar_type;

		inline vec() { }
		inline vec(eZero) { set_zero(); }

		explicit inline vec(T val) { set(val); }
		inline vec(T v0, T v1) { set(v0, v1); }
		inline vec(T v0, T v1, T v2) { set(v0, v1, v2); }
		inline vec(T v0, T v1, T v2, T v3) { set(v0, v1, v2, v3); }
		inline vec(const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] = other.m_v[i]; }
		template <uint32_t OtherN, typename OtherT> inline vec(const vec<OtherN, OtherT>& other) { set(other); }

		inline T operator[](uint32_t i) const { assert(i < N); return m_v[i]; }
		inline T& operator[](uint32_t i) { assert(i < N); return m_v[i]; }

		inline T getX() const { return m_v[0]; }
		inline T getY() const { static_assert(N >= 2, "N too small"); return m_v[1]; }
		inline T getZ() const { static_assert(N >= 3, "N too small"); return m_v[2]; }
		inline T getW() const { static_assert(N >= 4, "N too small"); return m_v[3]; }

		inline bool operator==(const vec& rhs) const { for (uint32_t i = 0; i < N; i++) if (m_v[i] != rhs.m_v[i]) return false;	return true; }
		inline bool operator<(const vec& rhs) const { for (uint32_t i = 0; i < N; i++) { if (m_v[i] < rhs.m_v[i]) return true; else if (m_v[i] != rhs.m_v[i]) return false; } return false; }

		inline void set_zero() { for (uint32_t i = 0; i < N; i++) m_v[i] = 0; }
		inline void clear() { set_zero(); }

		template <uint32_t OtherN, typename OtherT>
		inline vec& set(const vec<OtherN, OtherT>& other)
		{
			uint32_t i;
			if ((const void*)(&other) == (const void*)(this))
				return *this;
			const uint32_t m = minimum(OtherN, N);
			for (i = 0; i < m; i++)
				m_v[i] = static_cast<T>(other[i]);
			for (; i < N; i++)
				m_v[i] = 0;
			return *this;
		}

		inline vec& set_component(uint32_t index, T val) { assert(index < N); m_v[index] = val; return *this; }
		inline vec& set(T val) { for (uint32_t i = 0; i < N; i++) m_v[i] = val; return *this; }
		inline void clear_elements(uint32_t s, uint32_t e) { assert(e <= N); for (uint32_t i = s; i < e; i++) m_v[i] = 0; }

		inline vec& set(T v0, T v1)
		{
			m_v[0] = v0;
			if (N >= 2)
			{
				m_v[1] = v1;
				clear_elements(2, N);
			}
			return *this;
		}

		inline vec& set(T v0, T v1, T v2)
		{
			m_v[0] = v0;
			if (N >= 2)
			{
				m_v[1] = v1;
				if (N >= 3)
				{
					m_v[2] = v2;
					clear_elements(3, N);
				}
			}
			return *this;
		}

		inline vec& set(T v0, T v1, T v2, T v3)
		{
			m_v[0] = v0;
			if (N >= 2)
			{
				m_v[1] = v1;
				if (N >= 3)
				{
					m_v[2] = v2;

					if (N >= 4)
					{
						m_v[3] = v3;
						clear_elements(5, N);
					}
				}
			}
			return *this;
		}

		inline vec& operator=(const vec& rhs) { if (this != &rhs) for (uint32_t i = 0; i < N; i++) m_v[i] = rhs.m_v[i]; return *this; }
		template <uint32_t OtherN, typename OtherT> inline vec& operator=(const vec<OtherN, OtherT>& rhs) { set(rhs); return *this; }

		inline const T* get_ptr() const { return reinterpret_cast<const T*>(&m_v[0]); }
		inline T* get_ptr() { return reinterpret_cast<T*>(&m_v[0]); }

		inline vec operator- () const { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = -m_v[i]; return res; }
		inline vec operator+ () const { return *this; }
		inline vec& operator+= (const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] += other.m_v[i]; return *this; }
		inline vec& operator-= (const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] -= other.m_v[i]; return *this; }
		inline vec& operator/= (const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] /= other.m_v[i]; return *this; }
		inline vec& operator*=(const vec& other) { for (uint32_t i = 0; i < N; i++) m_v[i] *= other.m_v[i]; return *this; }
		inline vec& operator/= (T s) { for (uint32_t i = 0; i < N; i++) m_v[i] /= s; return *this; }
		inline vec& operator*= (T s) { for (uint32_t i = 0; i < N; i++) m_v[i] *= s; return *this; }

		friend inline vec operator+(const vec& lhs, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] + rhs.m_v[i]; return res; }
		friend inline vec operator-(const vec& lhs, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] - rhs.m_v[i]; return res; }
		friend inline vec operator*(const vec& lhs, T val) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] * val; return res; }
		friend inline vec operator*(T val, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = val * rhs.m_v[i]; return res; }
		friend inline vec operator/(const vec& lhs, T val) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] / val; return res; }
		friend inline vec operator/(const vec& lhs, const vec& rhs) { vec res; for (uint32_t i = 0; i < N; i++) res.m_v[i] = lhs.m_v[i] / rhs.m_v[i]; return res; }

		static inline T dot_product(const vec& lhs, const vec& rhs) { T res = lhs.m_v[0] * rhs.m_v[0]; for (uint32_t i = 1; i < N; i++) res += lhs.m_v[i] * rhs.m_v[i]; return res; }

		inline T dot(const vec& rhs) const { return dot_product(*this, rhs); }

		inline T norm() const { return dot_product(*this, *this); }
		inline T length() const { return sqrt(norm()); }

		inline T squared_distance(const vec& other) const { T d2 = 0; for (uint32_t i = 0; i < N; i++) { T d = m_v[i] - other.m_v[i]; d2 += d * d; } return d2; }
		inline double squared_distance_d(const vec& other) const { double d2 = 0; for (uint32_t i = 0; i < N; i++) { double d = (double)m_v[i] - (double)other.m_v[i]; d2 += d * d; } return d2; }

		inline T distance(const vec& other) const { return static_cast<T>(sqrt(squared_distance(other))); }
		inline double distance_d(const vec& other) const { return sqrt(squared_distance_d(other)); }

		inline vec& normalize_in_place() { T len = length(); if (len != 0.0f) *this *= (1.0f / len); return *this; }

		inline vec& clamp(T l, T h)
		{
			for (uint32_t i = 0; i < N; i++)
				m_v[i] = basisu::clamp(m_v[i], l, h);
			return *this;
		}

		static vec component_min(const vec& a, const vec& b)
		{
			vec res;
			for (uint32_t i = 0; i < N; i++)
				res[i] = minimum(a[i], b[i]);
			return res;
		}

		static vec component_max(const vec& a, const vec& b)
		{
			vec res;
			for (uint32_t i = 0; i < N; i++)
				res[i] = maximum(a[i], b[i]);
			return res;
		}
	};

	typedef vec<4, double> vec4D;
	typedef vec<3, double> vec3D;
	typedef vec<2, double> vec2D;
	typedef vec<1, double> vec1D;

	typedef vec<4, float> vec4F;
	typedef vec<3, float> vec3F;
	typedef vec<2, float> vec2F;
	typedef vec<1, float> vec1F;

	typedef vec<16, float> vec16F;

	class color_rgba
	{
	public:
		union
		{
			uint8_t m_comps[4];

			struct
			{
				uint8_t r;
				uint8_t g;
				uint8_t b;
				uint8_t a;
			};
		};

		inline color_rgba()
		{
			static_assert(sizeof(*this) == 4, "sizeof(*this) != 4");
		}

		inline color_rgba(int y)
		{
			set(y);
		}

		inline color_rgba(int y, int na)
		{
			set(y, na);
		}

		inline color_rgba(int sr, int sg, int sb, int sa)
		{
			set(sr, sg, sb, sa);
		}

		inline color_rgba(eNoClamp, int sr, int sg, int sb, int sa)
		{
			set_noclamp_rgba((uint8_t)sr, (uint8_t)sg, (uint8_t)sb, (uint8_t)sa);
		}

		inline color_rgba& set_noclamp_y(int y)
		{
			m_comps[0] = (uint8_t)y;
			m_comps[1] = (uint8_t)y;
			m_comps[2] = (uint8_t)y;
			m_comps[3] = (uint8_t)255;
			return *this;
		}

		inline color_rgba& set_noclamp_rgba(int sr, int sg, int sb, int sa)
		{
			m_comps[0] = (uint8_t)sr;
			m_comps[1] = (uint8_t)sg;
			m_comps[2] = (uint8_t)sb;
			m_comps[3] = (uint8_t)sa;
			return *this;
		}

		inline color_rgba& set(int y)
		{
			m_comps[0] = static_cast<uint8_t>(clamp<int>(y, 0, 255));
			m_comps[1] = m_comps[0];
			m_comps[2] = m_comps[0];
			m_comps[3] = 255;
			return *this;
		}

		inline color_rgba& set(int y, int na)
		{
			m_comps[0] = static_cast<uint8_t>(clamp<int>(y, 0, 255));
			m_comps[1] = m_comps[0];
			m_comps[2] = m_comps[0];
			m_comps[3] = static_cast<uint8_t>(clamp<int>(na, 0, 255));
			return *this;
		}

		inline color_rgba& set(int sr, int sg, int sb, int sa)
		{
			m_comps[0] = static_cast<uint8_t>(clamp<int>(sr, 0, 255));
			m_comps[1] = static_cast<uint8_t>(clamp<int>(sg, 0, 255));
			m_comps[2] = static_cast<uint8_t>(clamp<int>(sb, 0, 255));
			m_comps[3] = static_cast<uint8_t>(clamp<int>(sa, 0, 255));
			return *this;
		}

		inline color_rgba& set_rgb(int sr, int sg, int sb)
		{
			m_comps[0] = static_cast<uint8_t>(clamp<int>(sr, 0, 255));
			m_comps[1] = static_cast<uint8_t>(clamp<int>(sg, 0, 255));
			m_comps[2] = static_cast<uint8_t>(clamp<int>(sb, 0, 255));
			return *this;
		}

		inline color_rgba& set_rgb(const color_rgba& other)
		{
			r = other.r;
			g = other.g;
			b = other.b;
			return *this;
		}

		inline const uint8_t& operator[] (uint32_t index) const { assert(index < 4); return m_comps[index]; }
		inline uint8_t& operator[] (uint32_t index) { assert(index < 4); return m_comps[index]; }

		inline void clear()
		{
			m_comps[0] = 0;
			m_comps[1] = 0;
			m_comps[2] = 0;
			m_comps[3] = 0;
		}

		inline bool operator== (const color_rgba& rhs) const
		{
			if (m_comps[0] != rhs.m_comps[0]) return false;
			if (m_comps[1] != rhs.m_comps[1]) return false;
			if (m_comps[2] != rhs.m_comps[2]) return false;
			if (m_comps[3] != rhs.m_comps[3]) return false;
			return true;
		}

		inline bool operator!= (const color_rgba& rhs) const
		{
			return !(*this == rhs);
		}

		inline bool operator<(const color_rgba& rhs) const
		{
			for (int i = 0; i < 4; i++)
			{
				if (m_comps[i] < rhs.m_comps[i])
					return true;
				else if (m_comps[i] != rhs.m_comps[i])
					return false;
			}
			return false;
		}

		inline int get_601_luma() const { return (19595U * m_comps[0] + 38470U * m_comps[1] + 7471U * m_comps[2] + 32768U) >> 16U; }
		inline int get_709_luma() const { return (13938U * m_comps[0] + 46869U * m_comps[1] + 4729U * m_comps[2] + 32768U) >> 16U; }
		inline int get_luma(bool luma_601) const { return luma_601 ? get_601_luma() : get_709_luma(); }

		static color_rgba comp_min(const color_rgba& a, const color_rgba& b) { return color_rgba(basisu::minimum(a[0], b[0]), basisu::minimum(a[1], b[1]), basisu::minimum(a[2], b[2]), basisu::minimum(a[3], b[3])); }
		static color_rgba comp_max(const color_rgba& a, const color_rgba& b) { return color_rgba(basisu::maximum(a[0], b[0]), basisu::maximum(a[1], b[1]), basisu::maximum(a[2], b[2]), basisu::maximum(a[3], b[3])); }
	};

	typedef std::vector<color_rgba> color_rgba_vec;

	const color_rgba g_black_color(0, 0, 0, 255);
	const color_rgba g_black_trans_color(0, 0, 0, 0);
	const color_rgba g_white_color(255, 255, 255, 255);

	//--------------------------------------------------------------------------------------------------------------------------

	// Simple 32-bit 2D image class

	class image
	{
	public:
		image() :
			m_width(0), m_height(0), m_pitch(0)
		{
		}

		image(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX) :
			m_width(0), m_height(0), m_pitch(0)
		{
			resize(w, h, p);
		}

		image(const uint8_t* pImage, uint32_t width, uint32_t height, uint32_t comps) :
			m_width(0), m_height(0), m_pitch(0)
		{
			init(pImage, width, height, comps);
		}

		image(const image& other) :
			m_width(0), m_height(0), m_pitch(0)
		{
			*this = other;
		}

		image& swap(image& other)
		{
			std::swap(m_width, other.m_width);
			std::swap(m_height, other.m_height);
			std::swap(m_pitch, other.m_pitch);
			m_pixels.swap(other.m_pixels);
			return *this;
		}

		image& operator= (const image& rhs)
		{
			if (this != &rhs)
			{
				m_width = rhs.m_width;
				m_height = rhs.m_height;
				m_pitch = rhs.m_pitch;
				m_pixels = rhs.m_pixels;
			}
			return *this;
		}

		image& clear()
		{
			m_width = 0;
			m_height = 0;
			m_pitch = 0;
			clear_vector(m_pixels);
			return *this;
		}

		image& resize(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX, const color_rgba& background = g_black_color)
		{
			return crop(w, h, p, background);
		}

		image& set_all(const color_rgba& c)
		{
			for (uint32_t i = 0; i < m_pixels.size(); i++)
				m_pixels[i] = c;
			return *this;
		}

		void init(const uint8_t* pImage, uint32_t width, uint32_t height, uint32_t comps)
		{
			assert(comps >= 1 && comps <= 4);

			resize(width, height);

			for (uint32_t y = 0; y < height; y++)
			{
				for (uint32_t x = 0; x < width; x++)
				{
					const uint8_t* pSrc = &pImage[(x + y * width) * comps];
					color_rgba& dst = (*this)(x, y);

					if (comps == 1)
					{
						dst.r = pSrc[0];
						dst.g = pSrc[0];
						dst.b = pSrc[0];
						dst.a = 255;
					}
					else if (comps == 2)
					{
						dst.r = pSrc[0];
						dst.g = pSrc[0];
						dst.b = pSrc[0];
						dst.a = pSrc[1];
					}
					else
					{
						dst.r = pSrc[0];
						dst.g = pSrc[1];
						dst.b = pSrc[2];
						if (comps == 4)
							dst.a = pSrc[3];
						else
							dst.a = 255;
					}
				}
			}
		}

		image& fill_box(uint32_t x, uint32_t y, uint32_t w, uint32_t h, const color_rgba& c)
		{
			for (uint32_t iy = 0; iy < h; iy++)
				for (uint32_t ix = 0; ix < w; ix++)
					set_clipped(x + ix, y + iy, c);
			return *this;
		}

		image& fill_box_alpha(uint32_t x, uint32_t y, uint32_t w, uint32_t h, const color_rgba& c)
		{
			for (uint32_t iy = 0; iy < h; iy++)
				for (uint32_t ix = 0; ix < w; ix++)
					set_clipped_alpha(x + ix, y + iy, c);
			return *this;
		}

		image& crop_dup_borders(uint32_t w, uint32_t h)
		{
			const uint32_t orig_w = m_width, orig_h = m_height;

			crop(w, h);

			if (orig_w && orig_h)
			{
				if (m_width > orig_w)
				{
					for (uint32_t x = orig_w; x < m_width; x++)
						for (uint32_t y = 0; y < m_height; y++)
							set_clipped(x, y, get_clamped(minimum(x, orig_w - 1U), minimum(y, orig_h - 1U)));
				}

				if (m_height > orig_h)
				{
					for (uint32_t y = orig_h; y < m_height; y++)
						for (uint32_t x = 0; x < m_width; x++)
							set_clipped(x, y, get_clamped(minimum(x, orig_w - 1U), minimum(y, orig_h - 1U)));
				}
			}
			return *this;
		}

		image& crop(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX, const color_rgba& background = g_black_color, bool init_image = true)
		{
			if (p == UINT32_MAX)
				p = w;

			if ((w == m_width) && (m_height == h) && (m_pitch == p))
				return *this;

			if ((!w) || (!h) || (!p))
			{
				clear();
				return *this;
			}

			color_rgba_vec cur_state;
			cur_state.swap(m_pixels);

			m_pixels.resize(p * h);

			if (init_image)
			{
				if (m_width || m_height)
				{
					for (uint32_t y = 0; y < h; y++)
					{
						for (uint32_t x = 0; x < w; x++)
						{
							if ((x < m_width) && (y < m_height))
								m_pixels[x + y * p] = cur_state[x + y * m_pitch];
							else
								m_pixels[x + y * p] = background;
						}
					}
				}
				else
				{
					//m_pixels.set_all(background);
					for (uint32_t i = 0; i < m_pixels.size(); i++)
						m_pixels[i] = background;
				}
			}

			m_width = w;
			m_height = h;
			m_pitch = p;

			return *this;
		}

		inline const color_rgba& operator() (uint32_t x, uint32_t y) const { assert(x < m_width && y < m_height); return m_pixels[x + y * m_pitch]; }
		inline color_rgba& operator() (uint32_t x, uint32_t y) { assert(x < m_width && y < m_height); return m_pixels[x + y * m_pitch]; }

		inline const color_rgba& get_clamped(int x, int y) const { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }
		inline color_rgba& get_clamped(int x, int y) { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }

		inline const color_rgba& get_clamped_or_wrapped(int x, int y, bool wrap_u, bool wrap_v) const
		{
			x = wrap_u ? posmod(x, m_width) : clamp<int>(x, 0, m_width - 1);
			y = wrap_v ? posmod(y, m_height) : clamp<int>(y, 0, m_height - 1);
			return m_pixels[x + y * m_pitch];
		}

		inline color_rgba& get_clamped_or_wrapped(int x, int y, bool wrap_u, bool wrap_v)
		{
			x = wrap_u ? posmod(x, m_width) : clamp<int>(x, 0, m_width - 1);
			y = wrap_v ? posmod(y, m_height) : clamp<int>(y, 0, m_height - 1);
			return m_pixels[x + y * m_pitch];
		}

		inline image& set_clipped(int x, int y, const color_rgba& c)
		{
			if ((static_cast<uint32_t>(x) < m_width) && (static_cast<uint32_t>(y) < m_height))
				(*this)(x, y) = c;
			return *this;
		}

		inline image& set_clipped_alpha(int x, int y, const color_rgba& c)
		{
			if ((static_cast<uint32_t>(x) < m_width) && (static_cast<uint32_t>(y) < m_height))
				(*this)(x, y).m_comps[3] = c.m_comps[3];
			return *this;
		}

		// Very straightforward blit with full clipping. Not fast, but it works.
		image& blit(const image& src, int src_x, int src_y, int src_w, int src_h, int dst_x, int dst_y)
		{
			for (int y = 0; y < src_h; y++)
			{
				const int sy = src_y + y;
				if (sy < 0)
					continue;
				else if (sy >= (int)src.get_height())
					break;

				for (int x = 0; x < src_w; x++)
				{
					const int sx = src_x + x;
					if (sx < 0)
						continue;
					else if (sx >= (int)src.get_height())
						break;

					set_clipped(dst_x + x, dst_y + y, src(sx, sy));
				}
			}

			return *this;
		}

		const image& extract_block_clamped(color_rgba* pDst, uint32_t src_x, uint32_t src_y, uint32_t w, uint32_t h) const
		{
			if (((src_x + w) > m_width) || ((src_y + h) > m_height))
			{
				// Slower clamping case
				for (uint32_t y = 0; y < h; y++)
					for (uint32_t x = 0; x < w; x++)
						*pDst++ = get_clamped(src_x + x, src_y + y);
			}
			else
			{
				const color_rgba* pSrc = &m_pixels[src_x + src_y * m_pitch];

				for (uint32_t y = 0; y < h; y++)
				{
					memcpy(pDst, pSrc, w * sizeof(color_rgba));
					pSrc += m_pitch;
					pDst += w;
				}
			}

			return *this;
		}

		image& set_block_clipped(const color_rgba* pSrc, uint32_t dst_x, uint32_t dst_y, uint32_t w, uint32_t h)
		{
			for (uint32_t y = 0; y < h; y++)
				for (uint32_t x = 0; x < w; x++)
					set_clipped(dst_x + x, dst_y + y, *pSrc++);
			return *this;
		}

		inline uint32_t get_width() const { return m_width; }
		inline uint32_t get_height() const { return m_height; }
		inline uint32_t get_pitch() const { return m_pitch; }
		inline uint32_t get_total_pixels() const { return m_width * m_height; }

		inline uint32_t get_block_width(uint32_t w) const { return (m_width + (w - 1)) / w; }
		inline uint32_t get_block_height(uint32_t h) const { return (m_height + (h - 1)) / h; }
		inline uint32_t get_total_blocks(uint32_t w, uint32_t h) const { return get_block_width(w) * get_block_height(h); }

		inline const color_rgba_vec& get_pixels() const { return m_pixels; }
		inline color_rgba_vec& get_pixels() { return m_pixels; }

		inline const color_rgba* get_ptr() const { return &m_pixels[0]; }
		inline color_rgba* get_ptr() { return &m_pixels[0]; }

		bool has_alpha() const
		{
			for (uint32_t y = 0; y < m_height; ++y)
				for (uint32_t x = 0; x < m_width; ++x)
					if ((*this)(x, y).a < 255)
						return true;

			return false;
		}

		image& set_alpha(uint8_t a)
		{
			for (uint32_t y = 0; y < m_height; ++y)
				for (uint32_t x = 0; x < m_width; ++x)
					(*this)(x, y).a = a;
			return *this;
		}

		image& flip_y()
		{
			for (uint32_t y = 0; y < m_height / 2; ++y)
				for (uint32_t x = 0; x < m_width; ++x)
					std::swap((*this)(x, y), (*this)(x, m_height - 1 - y));
			return *this;
		}

		// TODO: There are many ways to do this, not sure this is the best way.
		image& renormalize_normal_map()
		{
			for (uint32_t y = 0; y < m_height; y++)
			{
				for (uint32_t x = 0; x < m_width; x++)
				{
					color_rgba& c = (*this)(x, y);
					if ((c.r == 128) && (c.g == 128) && (c.b == 128))
						continue;

					vec3F v(c.r, c.g, c.b);
					v = (v * (2.0f / 255.0f)) - vec3F(1.0f);
					v.clamp(-1.0f, 1.0f);

					float length = v.length();
					const float cValidThresh = .077f;
					if (length < cValidThresh)
					{
						c.set(128, 128, 128, c.a);
					}
					else if (fabs(length - 1.0f) > cValidThresh)
					{
						if (length)
							v /= length;

						for (uint32_t i = 0; i < 3; i++)
							c[i] = static_cast<uint8_t>(clamp<float>(floorf((v[i] + 1.0f) * 255.0f * .5f + .5f), 0.0f, 255.0f));

						if ((c.g == 128) && (c.r == 128))
						{
							if (c.b < 128)
								c.b = 0;
							else
								c.b = 255;
						}
					}
				}
			}
			return *this;
		}

	private:
		uint32_t m_width, m_height, m_pitch;  // all in pixels
		color_rgba_vec m_pixels;
	};

	// Float images

	typedef std::vector<vec4F> vec4F_vec;

	class imagef
	{
	public:
		imagef() :
			m_width(0), m_height(0), m_pitch(0)
		{
		}

		imagef(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX) :
			m_width(0), m_height(0), m_pitch(0)
		{
			resize(w, h, p);
		}

		imagef(const imagef& other) :
			m_width(0), m_height(0), m_pitch(0)
		{
			*this = other;
		}

		imagef& swap(imagef& other)
		{
			std::swap(m_width, other.m_width);
			std::swap(m_height, other.m_height);
			std::swap(m_pitch, other.m_pitch);
			m_pixels.swap(other.m_pixels);
			return *this;
		}

		imagef& operator= (const imagef& rhs)
		{
			if (this != &rhs)
			{
				m_width = rhs.m_width;
				m_height = rhs.m_height;
				m_pitch = rhs.m_pitch;
				m_pixels = rhs.m_pixels;
			}
			return *this;
		}

		imagef& clear()
		{
			m_width = 0;
			m_height = 0;
			m_pitch = 0;
			clear_vector(m_pixels);
			return *this;
		}

		imagef& set(const image& src, const vec4F& scale = vec4F(1), const vec4F& bias = vec4F(0))
		{
			const uint32_t width = src.get_width();
			const uint32_t height = src.get_height();

			resize(width, height);

			for (int y = 0; y < (int)height; y++)
			{
				for (uint32_t x = 0; x < width; x++)
				{
					const color_rgba& src_pixel = src(x, y);
					(*this)(x, y).set((float)src_pixel.r * scale[0] + bias[0], (float)src_pixel.g * scale[1] + bias[1], (float)src_pixel.b * scale[2] + bias[2], (float)src_pixel.a * scale[3] + bias[3]);
				}
			}

			return *this;
		}

		imagef& resize(const imagef& other, uint32_t p = UINT32_MAX, const vec4F& background = vec4F(0, 0, 0, 1))
		{
			return resize(other.get_width(), other.get_height(), p, background);
		}

		imagef& resize(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX, const vec4F& background = vec4F(0, 0, 0, 1))
		{
			return crop(w, h, p, background);
		}

		imagef& set_all(const vec4F& c)
		{
			for (uint32_t i = 0; i < m_pixels.size(); i++)
				m_pixels[i] = c;
			return *this;
		}

		imagef& fill_box(uint32_t x, uint32_t y, uint32_t w, uint32_t h, const vec4F& c)
		{
			for (uint32_t iy = 0; iy < h; iy++)
				for (uint32_t ix = 0; ix < w; ix++)
					set_clipped(x + ix, y + iy, c);
			return *this;
		}

		imagef& crop(uint32_t w, uint32_t h, uint32_t p = UINT32_MAX, const vec4F& background = vec4F(0, 0, 0, 1))
		{
			if (p == UINT32_MAX)
				p = w;

			if ((w == m_width) && (m_height == h) && (m_pitch == p))
				return *this;

			if ((!w) || (!h) || (!p))
			{
				clear();
				return *this;
			}

			vec4F_vec cur_state;
			cur_state.swap(m_pixels);

			m_pixels.resize(p * h);

			for (uint32_t y = 0; y < h; y++)
			{
				for (uint32_t x = 0; x < w; x++)
				{
					if ((x < m_width) && (y < m_height))
						m_pixels[x + y * p] = cur_state[x + y * m_pitch];
					else
						m_pixels[x + y * p] = background;
				}
			}

			m_width = w;
			m_height = h;
			m_pitch = p;

			return *this;
		}

		inline const vec4F& operator() (uint32_t x, uint32_t y) const { assert(x < m_width && y < m_height); return m_pixels[x + y * m_pitch]; }
		inline vec4F& operator() (uint32_t x, uint32_t y) { assert(x < m_width && y < m_height); return m_pixels[x + y * m_pitch]; }

		inline const vec4F& get_clamped(int x, int y) const { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }
		inline vec4F& get_clamped(int x, int y) { return (*this)(clamp<int>(x, 0, m_width - 1), clamp<int>(y, 0, m_height - 1)); }

		inline const vec4F& get_clamped_or_wrapped(int x, int y, bool wrap_u, bool wrap_v) const
		{
			x = wrap_u ? posmod(x, m_width) : clamp<int>(x, 0, m_width - 1);
			y = wrap_v ? posmod(y, m_height) : clamp<int>(y, 0, m_height - 1);
			return m_pixels[x + y * m_pitch];
		}

		inline vec4F& get_clamped_or_wrapped(int x, int y, bool wrap_u, bool wrap_v)
		{
			x = wrap_u ? posmod(x, m_width) : clamp<int>(x, 0, m_width - 1);
			y = wrap_v ? posmod(y, m_height) : clamp<int>(y, 0, m_height - 1);
			return m_pixels[x + y * m_pitch];
		}

		inline imagef& set_clipped(int x, int y, const vec4F& c)
		{
			if ((static_cast<uint32_t>(x) < m_width) && (static_cast<uint32_t>(y) < m_height))
				(*this)(x, y) = c;
			return *this;
		}

		// Very straightforward blit with full clipping. Not fast, but it works.
		imagef& blit(const imagef& src, int src_x, int src_y, int src_w, int src_h, int dst_x, int dst_y)
		{
			for (int y = 0; y < src_h; y++)
			{
				const int sy = src_y + y;
				if (sy < 0)
					continue;
				else if (sy >= (int)src.get_height())
					break;

				for (int x = 0; x < src_w; x++)
				{
					const int sx = src_x + x;
					if (sx < 0)
						continue;
					else if (sx >= (int)src.get_height())
						break;

					set_clipped(dst_x + x, dst_y + y, src(sx, sy));
				}
			}

			return *this;
		}

		const imagef& extract_block_clamped(vec4F* pDst, uint32_t src_x, uint32_t src_y, uint32_t w, uint32_t h) const
		{
			for (uint32_t y = 0; y < h; y++)
				for (uint32_t x = 0; x < w; x++)
					*pDst++ = get_clamped(src_x + x, src_y + y);
			return *this;
		}

		imagef& set_block_clipped(const vec4F* pSrc, uint32_t dst_x, uint32_t dst_y, uint32_t w, uint32_t h)
		{
			for (uint32_t y = 0; y < h; y++)
				for (uint32_t x = 0; x < w; x++)
					set_clipped(dst_x + x, dst_y + y, *pSrc++);
			return *this;
		}

		inline uint32_t get_width() const { return m_width; }
		inline uint32_t get_height() const { return m_height; }
		inline uint32_t get_pitch() const { return m_pitch; }
		inline uint32_t get_total_pixels() const { return m_width * m_height; }

		inline uint32_t get_block_width(uint32_t w) const { return (m_width + (w - 1)) / w; }
		inline uint32_t get_block_height(uint32_t h) const { return (m_height + (h - 1)) / h; }
		inline uint32_t get_total_blocks(uint32_t w, uint32_t h) const { return get_block_width(w) * get_block_height(h); }

		inline const vec4F_vec& get_pixels() const { return m_pixels; }
		inline vec4F_vec& get_pixels() { return m_pixels; }

		inline const vec4F* get_ptr() const { return &m_pixels[0]; }
		inline vec4F* get_ptr() { return &m_pixels[0]; }

		bool clean_pixels()
		{
			bool status = true;

			for (uint32_t y = 0; y < m_height; y++)
			{
				for (uint32_t x = 0; x < m_width; x++)
				{
					vec4F& c = (*this)(x, y);

					for (uint32_t s = 0; s < 4; s++)
					{
						float& p = c[s];

						if ((std::isnan(p)) || (std::isinf(p)))
						{
							p = 0.0f;
							status = false;
						}
						else
						{
							const float o = p;
														
							p = basisu::minimumf(MAX_HALF_FLOAT, p);

							if (p != o)
								status = false;
						}
					}
				}
			}

			return status;
		}

	private:
		uint32_t m_width, m_height, m_pitch;  // all in pixels
		vec4F_vec m_pixels;
	};

	enum
	{
		cImageSaveGrayscale = 1,
		cImageSaveIgnoreAlpha = 2
	};

	bool load_png(const uint8_t* pBuf, size_t buf_size, image& img, const char* pFilename)
	{
		(void)pFilename;

		if (!buf_size)
			return false;

		unsigned err = 0, w = 0, h = 0;

		std::vector<uint8_t> out;
		err = lodepng::decode(out, w, h, pBuf, buf_size);
		if ((err != 0) || (!w) || (!h))
			return false;

		if (out.size() != (w * h * 4))
			return false;

		img.resize(w, h);

		memcpy(img.get_ptr(), &out[0], out.size());

		return true;
	}

	bool save_png(const char* pFilename, const image& img, uint32_t image_save_flags, uint32_t grayscale_comp)
	{
		if (!img.get_total_pixels())
			return false;

		std::vector<uint8_t> out;
		unsigned err = 0;

		if (image_save_flags & cImageSaveGrayscale)
		{
			uint8_vec g_pixels(img.get_width() * img.get_height());
			uint8_t* pDst = &g_pixels[0];

			for (uint32_t y = 0; y < img.get_height(); y++)
				for (uint32_t x = 0; x < img.get_width(); x++)
					*pDst++ = img(x, y)[grayscale_comp];

			err = lodepng::encode(out, (const uint8_t*)&g_pixels[0], img.get_width(), img.get_height(), LCT_GREY, 8);
		}
		else
		{
			bool has_alpha = img.has_alpha();
			if ((!has_alpha) || ((image_save_flags & cImageSaveIgnoreAlpha) != 0))
			{
				const uint64_t total_bytes = (uint64_t)img.get_width() * 3U * (uint64_t)img.get_height();
				if (total_bytes > INT_MAX)
					return false;
				uint8_vec rgb_pixels(static_cast<size_t>(total_bytes));
				uint8_t* pDst = &rgb_pixels[0];

				for (uint32_t y = 0; y < img.get_height(); y++)
				{
					for (uint32_t x = 0; x < img.get_width(); x++)
					{
						const color_rgba& c = img(x, y);
						pDst[0] = c.r;
						pDst[1] = c.g;
						pDst[2] = c.b;
						pDst += 3;
					}
				}

				err = lodepng::encode(out, (const uint8_t*)&rgb_pixels[0], img.get_width(), img.get_height(), LCT_RGB, 8);
			}
			else
			{
				err = lodepng::encode(out, (const uint8_t*)img.get_ptr(), img.get_width(), img.get_height(), LCT_RGBA, 8);
			}
		}

		err = lodepng::save_file(out, std::string(pFilename));
		if (err)
			return false;

		return true;
	}

#pragma pack(push, 1)
	struct png_hdr_chunk
	{
		uint32_t m_png_len;		// size of chunk data (only)
		uint32_t m_png_type;	// chunk ID

		// start of chunk data
		uint8_t m_shift_amount;
		uint8_t m_remap_table[256];
		// end of chunk data

		uint32_t m_png_crc;

		void init()
		{
			uint32_t len = (uint32_t)(sizeof(png_hdr_chunk) - sizeof(uint32_t) * 3);
			m_png_len = byteswap32(len);
			m_png_type = 'h' | ('d' << 8) | ('R' << 16) | ('a' << 24);
			m_png_crc = byteswap32(buminiz::mz_crc32(MZ_CRC32_INIT, (uint8_t*)&m_png_type, len + sizeof(uint32_t)));
		}

		bool check() const
		{
			uint32_t len = byteswap32(m_png_len);
			if (len != (sizeof(png_hdr_chunk) - sizeof(uint32_t) * 3))
				return false;
			if (m_png_type != ('h' | ('d' << 8) | ('R' << 16) | ('a' << 24)))
				return false;
			if (byteswap32(m_png_crc) != buminiz::mz_crc32(MZ_CRC32_INIT, (uint8_t*)&m_png_type, len + sizeof(uint32_t)))
				return false;
			return true;
		}
	};
#pragma pack(pop)

	bool load_rgb16_png(const char* pFilename, uint32_t& width, uint32_t& height, uint16_vec& img, png_hdr_chunk& hdr_chunk)
	{
		unsigned err = 0, w = 0, h = 0;

		clear_obj(hdr_chunk);

		uint8_vec in;
		if (lodepng::load_file(in, std::string(pFilename)) != 0)
			return false;

		if (!in.size())
			return false;

		std::vector<uint8_t> out;
		err = lodepng::decode(out, w, h, in, LCT_RGB, 16);
		if ((err != 0) || (!w) || (!h))
			return false;

		if (out.size() != (w * h * sizeof(uint16_t) * 3))
			return false;

		// Find hdRa chunk
		std::vector<std::string> names[3];
		std::vector<std::vector<unsigned char> > chunks[3];

		err = lodepng::getChunks(names, chunks, in);
		if (err != 0)
			return false;

		bool hdr_chunk_found = false;
		for (uint32_t i = 0; i < names[0].size(); i++)
		{
			if (names[0][i] == "hdRa")
			{
				if (chunks[0][i].size() == sizeof(png_hdr_chunk))
				{
					memcpy(&hdr_chunk, &chunks[0][i][0], sizeof(png_hdr_chunk));
					if (!hdr_chunk.check())
						return false;

					hdr_chunk_found = true;
				}
			}
		}

		width = w;
		height = h;

		img.resize(width * height * 3);
		for (uint32_t y = 0; y < height; y++)
		{
			for (uint32_t x = 0; x < width; x++)
			{
				const uint8_t* pPixel = &out[(x + y * width) * 6];

				img[(x + y * width) * 3 + 0] = (uint16_t)((pPixel[0] << 8) | pPixel[1]);
				img[(x + y * width) * 3 + 1] = (uint16_t)((pPixel[2] << 8) | pPixel[3]);
				img[(x + y * width) * 3 + 2] = (uint16_t)((pPixel[4] << 8) | pPixel[5]);
			}
		}

		return true;
	}

	bool save_rgb16_png(const char* pFilename, uint32_t width, uint32_t height, const uint16_vec& img, const png_hdr_chunk& hdr_chunk)
	{
		const uint32_t total_pixels = width * height;
		if (!total_pixels)
		{
			assert(0);
			return false;
		}

		if (img.size() != total_pixels * 3)
		{
			assert(0);
			return false;
		}

		std::vector<uint8_t> out;

		uint16_vec temp_img(img.size());
		for (uint32_t i = 0; i < temp_img.size(); i++)
		{
			uint32_t v = img[i];
			temp_img[i] = (uint16_t)((v >> 8) | (v << 8));
		}

		unsigned err = lodepng::encode(out, (const uint8_t*)&temp_img[0], width, height, LCT_RGB, 16);
		if (err)
			return false;

		// Add our private ancillary hdRa chunk to the PNG file data
		std::vector<std::vector<unsigned char> > chunks[3];
		chunks[0].push_back(std::vector<unsigned char>((unsigned char*)&hdr_chunk, (unsigned char*)&hdr_chunk + sizeof(hdr_chunk)));

		err = lodepng::insertChunks(out, chunks);
		if (err)
			return false;

		err = lodepng::save_file(out, std::string(pFilename));
		if (err)
			return false;

		return true;
	}

	enum
	{
		WRITE_EXR_LINEAR_HINT = 1, // hint for lossy comp. methods: exr_perceptual_treatment_t, logarithmic or linear, defaults to logarithmic
		WRITE_EXR_STORE_FLOATS = 2, // use 32-bit floats, otherwise it uses half floats
		WRITE_EXR_NO_COMPRESSION = 4 // no compression, otherwise it uses ZIP compression (16 scanlines per block)
	};

	bool read_exr(const char* pFilename, imagef& img, int& n_chans)
	{
		n_chans = 0;

		int width = 0, height = 0;
		float* out_rgba = nullptr;
		const char* err = nullptr;

		int status = LoadEXRWithLayer(&out_rgba, &width, &height, pFilename, nullptr, &err, &n_chans);
		if (status != 0)
		{
			fprintf(stderr, "Failed loading .EXR image \"%s\"! (TinyEXR error: %s)\n", pFilename, err);

			return false;
		}

		const uint32_t MAX_SUPPORTED_DIM = 65536;
		if ((width < 1) || (height < 1) || (width > (int)MAX_SUPPORTED_DIM) || (height > (int)MAX_SUPPORTED_DIM))
		{
			fprintf(stderr, "Invalid dimensions of .EXR image \"%s\"!\n", pFilename);

			free(out_rgba);
			return false;
		}

		img.resize(width, height);

		if (n_chans == 1)
		{
			const float* pSrc = out_rgba;
			vec4F* pDst = img.get_ptr();

			for (int y = 0; y < height; y++)
			{
				for (int x = 0; x < width; x++)
				{
					(*pDst)[0] = pSrc[0];
					(*pDst)[1] = pSrc[1];
					(*pDst)[2] = pSrc[2];
					(*pDst)[3] = 1.0f;

					pSrc += 4;
					++pDst;
				}
			}
		}
		else
		{
			memcpy(img.get_ptr(), out_rgba, sizeof(float) * 4 * img.get_total_pixels());
		}

		free(out_rgba);
		return true;
	}

	bool write_exr(const char* pFilename, imagef& img, uint32_t n_chans, uint32_t flags)
	{
		assert((n_chans == 1) || (n_chans == 3) || (n_chans == 4));

		const bool linear_hint = (flags & WRITE_EXR_LINEAR_HINT) != 0,
			store_float = (flags & WRITE_EXR_STORE_FLOATS) != 0,
			no_compression = (flags & WRITE_EXR_NO_COMPRESSION) != 0;

		const uint32_t width = img.get_width(), height = img.get_height();
		assert(width && height);

		if (!width || !height)
			return false;

		float_vec layers[4];
		float* image_ptrs[4];
		for (uint32_t c = 0; c < n_chans; c++)
		{
			layers[c].resize(width * height);
			image_ptrs[c] = &layers[c][0];// .get_ptr();
		}

		// ABGR
		int chan_order[4] = { 3, 2, 1, 0 };

		if (n_chans == 1)
		{
			// Y
			chan_order[0] = 0;
		}
		else if (n_chans == 3)
		{
			// BGR
			chan_order[0] = 2;
			chan_order[1] = 1;
			chan_order[2] = 0;
		}
		else if (n_chans != 4)
		{
			assert(0);
			return false;
		}

		for (uint32_t y = 0; y < height; y++)
		{
			for (uint32_t x = 0; x < width; x++)
			{
				const vec4F& p = img(x, y);

				for (uint32_t c = 0; c < n_chans; c++)
					layers[c][x + y * width] = p[chan_order[c]];
			} // x
		} // y

		EXRHeader header;
		InitEXRHeader(&header);

		EXRImage image;
		InitEXRImage(&image);

		image.num_channels = n_chans;
		image.images = (unsigned char**)image_ptrs;
		image.width = width;
		image.height = height;

		header.num_channels = n_chans;

		header.channels = (EXRChannelInfo*)calloc(header.num_channels, sizeof(EXRChannelInfo));

		// Must be (A)BGR order, since most of EXR viewers expect this channel order.
		for (uint32_t i = 0; i < n_chans; i++)
		{
			char c = 'Y';
			if (n_chans == 3)
				c = "BGR"[i];
			else if (n_chans == 4)
				c = "ABGR"[i];

			header.channels[i].name[0] = c;
			header.channels[i].name[1] = '\0';

			header.channels[i].p_linear = linear_hint;
		}

		header.pixel_types = (int*)calloc(header.num_channels, sizeof(int));
		header.requested_pixel_types = (int*)calloc(header.num_channels, sizeof(int));

		if (!no_compression)
			header.compression_type = TINYEXR_COMPRESSIONTYPE_ZIP;

		for (int i = 0; i < header.num_channels; i++)
		{
			// pixel type of input image
			header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;

			// pixel type of output image to be stored in .EXR
			header.requested_pixel_types[i] = store_float ? TINYEXR_PIXELTYPE_FLOAT : TINYEXR_PIXELTYPE_HALF;
		}

		const char* pErr_msg = nullptr;

		int ret = SaveEXRImageToFile(&image, &header, pFilename, &pErr_msg);
		if (ret != TINYEXR_SUCCESS)
		{
			fprintf(stderr, "Save EXR err: %s\n", pErr_msg);
			FreeEXRErrorMessage(pErr_msg);
		}

		free(header.channels);
		free(header.pixel_types);
		free(header.requested_pixel_types);

		return (ret == TINYEXR_SUCCESS);
	}

} // namespace basisu
