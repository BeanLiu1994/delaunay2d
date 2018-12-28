#pragma once

#include <string>
#include <vector>
#include <utility>
#include <numeric>
#include <array>
#define GenMatlabCmd
#define DELAUNAY_2D_USE_EIGEN

#ifdef DELAUNAY_2D_USE_EIGEN
#include <Eigen/Eigen>
#endif

namespace Delaunay2D
{
	typedef int Index;
	typedef double scalar;

	typedef std::pair<Index, Index> Edge;
	typedef std::array<scalar, 2> Point;
	typedef std::array<Index, 3> Triangle;

	inline bool operator<(const Point& p1, const Point& p2)
	{
		return p1[0] < p2[0] || (!(p2[0] < p1[0]) && p1[1] < p2[1]);
	}

	class Delaunay2D
	{
	private:
		void assign(const std::vector<Point>& Pts);
		void ForceSameOrientation();
		std::vector<Triangle> Tri;
	public:
		Delaunay2D(const std::vector<Point>& PtsToInsert);
		const std::vector<Triangle>& GetTri() const;

#ifdef DELAUNAY_2D_USE_EIGEN
	private:
		void assign(const Eigen::MatrixX2d& PtsToInsert);
	public:
		Delaunay2D(const Eigen::MatrixX2d& PtsToInsert);
		Eigen::MatrixX3i GetTri_Eigen() const;
#endif

#ifdef GenMatlabCmd
	private:
		std::string MatlabCmd;
	public:
		std::string GetMatlabCmd() const;
#endif
	};


	enum CircleCmpStatus
	{
		StrongReject,
		WeakReject,
		WeakAccept,
		Borderline
	};


	inline scalar dot(const Point& pts1, const Point& pts2)
	{
		return pts1[0] * pts2[0] + pts1[1] * pts2[1];
	}

	inline scalar dot(const Point& pts)
	{
		return dot(pts, pts);
	}

	inline Point operator-(const Point& pts1, const Point& pts2)
	{
		return Point{ pts1[0] - pts2[0], pts1[1] - pts2[1] };
	}

	/**
	 * from https://github.com/Bl4ckb0ne/delaunay-triangulation/blob/master/numeric.h
	 * @brief use of machine epsilon to compare floating-point values for equality
	 * http://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
	 */
	template<typename T>
	inline typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
		almost_equal(T x, T y, int ulp = 2)
	{
		// the machine epsilon has to be scaled to the magnitude of the values used
		// and multiplied by the desired precision in ULPs (units in the last place)
		return std::abs(x - y) <= std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
			// unless the result is subnormal
			|| std::abs(x - y) < std::numeric_limits<T>::min();
	}

	//compared by a sorted points (by x dim), here we treat x axis specially.
	//由于本作按照x排好序了，比较的时候会额外比较x的情况
	inline CircleCmpStatus isInsideTriangleCircumCircle(const Point& pts_now,
		const Point& p0, const Point& p1, const Point& p2)
	{
		// 解方程
		double p0_v = dot(p0);
		Point delta_p_val{ p0_v - dot(p1), p0_v - dot(p2) };

		Point
			delta_vec_p01 = (p0 - p1),
			delta_vec_p02 = (p0 - p2);

		double checker = 0.5 / (delta_vec_p01[0] * delta_vec_p02[1] - delta_vec_p01[1] * delta_vec_p02[0]);


		// 三点共线检查
		if (std::isinf(checker))
		{
			return WeakReject;
		}

		Point center = {
				checker * (delta_p_val[0] * delta_vec_p02[1] - delta_p_val[1] * delta_vec_p01[1]),
				checker * (delta_p_val[1] * delta_vec_p01[0] - delta_p_val[0] * delta_vec_p02[0])
		};

		double r = dot(p1 - center);
		

		// center 与 pts 的 x 比较, 假定已经按照x轴排序了
		// compare pts and center on x axis. assume that points are sorted by x axis.

		if (std::pow(pts_now[0] - center[0], 2) > r)
			return StrongReject;

		double len = dot(pts_now - center);

		if (almost_equal(len, r))
			return Borderline;

		if (len > r)
			return WeakReject;

		return WeakAccept;
	}


	// returns if no duplicate element is removed.
	template<typename ElemTy, typename Pred = std::greater<ElemTy>>
	bool RemoveDuplicateAndSort(std::vector<ElemTy>& ToProcess, Pred pred = Pred())
	{
		size_t n = ToProcess.size();
		std::sort(ToProcess.begin(), ToProcess.end(), pred);
		auto last = std::unique(ToProcess.begin(), ToProcess.end());
		ToProcess.erase(last, ToProcess.end());
		return (last - ToProcess.begin()) == n;
	}

	template<typename ElemTy, typename Pred = std::greater<ElemTy>>
	bool RemoveDuplicateAndSort(std::vector<ElemTy>& ToProcess, std::vector<int>& index, Pred pred = Pred())
	{
		size_t n = ToProcess.size();
		index.resize(ToProcess.size());
		std::iota(index.begin(), index.end(), 0);

		std::vector<std::pair<ElemTy, int>> combined;
		combined.reserve(ToProcess.size());
		for (int i = 0; i < ToProcess.size(); ++i)
			combined.emplace_back(ToProcess[i], index[i]);

		std::sort(combined.begin(), combined.end(),
			[&pred](std::pair<ElemTy, int>& v1, std::pair<ElemTy, int>& v2) {return pred(v1.first, v2.first); }
		);

		long long finish = std::unique(combined.begin(), combined.end(),
			[](std::pair<ElemTy, int>& v1, std::pair<ElemTy, int>& v2) {return v1.first == v2.first; }
		) - combined.begin();

		ToProcess.erase(ToProcess.begin() + finish, ToProcess.end());
		index.erase(index.begin() + finish, index.end());
		for (int i = 0; i < finish; ++i)
		{
			ToProcess[i] = combined[i].first;
			index[i] = combined[i].second;
		}
		return (finish) == n;
	}

#ifdef DELAUNAY_2D_USE_EIGEN
	template<typename scalar>
	std::vector<Point> MatToStdVec(const Eigen::Matrix<scalar, -1, 2>& In)
	{
		std::vector<Point> ret;
		ret.reserve(In.rows());
		for (int i = 0; i < In.rows(); ++i)
		{
			ret.emplace_back(Point{ In(i, 0), In(i, 1) });
		}
		return ret;
	}

#endif

}
