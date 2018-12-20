#pragma once

#include <Eigen/Eigen>
#include <string>
#include <utility>
//#define GenMatlabCmd

typedef std::pair<int, int> Edge;

class Delaunay2D
{
#ifdef GenMatlabCmd
	std::string MatlabCmd;
#endif
	void assign(const Eigen::MatrixX2d& PtsToInsert);
	void SameOrientation();
	Eigen::MatrixX3i Tri;
public:
	Delaunay2D(const Eigen::MatrixX2d& PtsToInsert);
	Eigen::MatrixX3i GetTri() const;
#ifdef GenMatlabCmd
	std::string GetMatlabCmd() const;
#endif
};


enum CircleCmpStatus
{
	StrongReject,
	WeakReject,
	WeakAccept
};

// compared by a sorted points (by x dim), here we treat x axis specially.
// 由于本作按照x排好序了，比较的时候会额外比较x的情况
inline CircleCmpStatus isInsideTriangleCircumCircle(const Eigen::RowVector2d& pts_now,
	const Eigen::RowVector2d& p0, const Eigen::RowVector2d& p1, const Eigen::RowVector2d& p2)
{
	// 解方程
	double p0_v = p0.dot(p0);
	Eigen::Vector2d delta_p_val(p0_v - p1.dot(p1), p0_v - p2.dot(p2));

	Eigen::RowVector2d 
		delta_vec_p01 = (p0 - p1).transpose(),
		delta_vec_p02 = (p0 - p2).transpose();

	double checker = 0.5 / (delta_vec_p01(0) * delta_vec_p02(1) - delta_vec_p01(1) * delta_vec_p02(0));

	Eigen::RowVector2d center(
		checker * (delta_p_val(0) * delta_vec_p02(1) - delta_p_val(1) * delta_vec_p01(1)),
		checker * (delta_p_val(1) * delta_vec_p01(0) - delta_p_val(0) * delta_vec_p02(0)));

	double r = (p1 - center).norm();

	// center 与 pts 的 x 比较
	// compare pts and center on x axis.
	if (pts_now(0) - center(0) > r)
		return StrongReject;

	if ((pts_now - center).norm() > r)
		return WeakReject;
	
	return WeakAccept;
}


// In --> sort
// dim: row  sort while treating each row as one part,  col sort while treating each col...
// referer: each row (or col) is expressed by one element inside. referer is the index of that value.
// pred: i didn't test this. 我没测试这个
//
enum dimension { row, col };
template<typename Derived, typename DerivedX, typename DerivedIX, typename Pred = std::less<typename Derived::Scalar>>
void sort(const Eigen::DenseBase<Derived>& In, const dimension dim, const int referer, Eigen::PlainObjectBase<DerivedX>& sorted, Eigen::PlainObjectBase<DerivedIX>& IX, Pred pred = Pred())
{
	static_assert(std::is_same<typename DerivedX::Scalar, typename Derived::Scalar>::value, "output should have same type as input.");
	typedef typename Derived::Scalar scalar;

	typename Eigen::DenseBase<Derived>::Index n;
	std::function<scalar(int)> access_to_sort_elem;
	if (dim == row)
	{
		n = In.rows();
		access_to_sort_elem = [&In, &referer](int ith_elem) {return In(ith_elem, referer); };
	}
	else if (dim == col)
	{
		n = In.cols();
		access_to_sort_elem = [&In, &referer](int ith_elem) {return In(referer, ith_elem); };
	}
	else
	{
		throw std::runtime_error("invalid option: dim.");
	}

	std::vector<std::pair<scalar, int>> ToSort(n);
	for (int i = 0; i < n; ++i)
	{
		ToSort[i].first = access_to_sort_elem(i);
		ToSort[i].second = i;
	}
	std::sort(ToSort.begin(), ToSort.end(),
		[&pred](const std::pair<scalar, int>& v1, const std::pair<scalar, int>& v2) {return pred(v1.first, v2.first); }
	);

	sorted.resizeLike(In);
	IX.resize(n);
	for (int i = 0; i < n; ++i)
	{
		IX(i) = ToSort[i].second;
		if (dim == row)
		{
			sorted.row(i) = In.row(ToSort[i].second);
		}
		else
		{
			sorted.col(i) = In.col(ToSort[i].second);
		}
	}
}
