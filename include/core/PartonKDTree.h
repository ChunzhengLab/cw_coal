#ifndef PARTON_KDTREE_H
#define PARTON_KDTREE_H

#include <cmath>
#include <nanoflann.hpp>
#include <vector>

#include "core/Particle.h"

class PartonCloud {
 public:
  PartonCloud(const std::vector<Parton*>& points) : m_points(points) {
  }

  inline size_t kdtree_get_point_count() const {
    return m_points.size();
  }

  inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
    if (dim == 0) return m_points[idx]->X();
    if (dim == 1) return m_points[idx]->Y();
    return m_points[idx]->Z();
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX&) const {
    return false;
  }

  Parton* get(size_t idx) const {
    return m_points[idx];
  }

 private:
  const std::vector<Parton*>& m_points;
};

using KDTree =
    nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PartonCloud>, PartonCloud, 3 /* dim */>;

class PartonKDTree {
 public:
  PartonKDTree(const std::vector<Parton*>& partons)
      : m_cloud(partons), m_index(3, m_cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10)) {
    m_index.buildIndex();
  }

  std::vector<std::pair<Parton*, double>> kNearestSearch(const Parton* query, size_t maxResults = 50) const {
    std::vector<std::pair<Parton*, double>> results;
    double queryPoint[3] = {query->X(), query->Y(), query->Z()};

    std::vector<size_t> ret_index(maxResults);
    std::vector<double> out_dist_sqr(maxResults);

    nanoflann::KNNResultSet<double> resultSet(maxResults);
    resultSet.init(ret_index.data(), out_dist_sqr.data());

    nanoflann::SearchParameters params;
    params.sorted = true;
    m_index.findNeighbors(resultSet, queryPoint, params);

    size_t found = resultSet.size();
    for (size_t i = 0; i < found; ++i) {
      Parton* p = m_cloud.get(ret_index[i]);
      if (p->IsUsed()) continue;
      results.emplace_back(p, std::sqrt(out_dist_sqr[i]));
    }

    return results;
  }

  std::vector<std::pair<Parton*, double>> FindNeighbors(const Parton* query, size_t maxResults = 50) const {
    return kNearestSearch(query, maxResults);
  }

  Parton* FindNearestOpposite(const Parton* query) const {
    double queryPoint[3] = {query->X(), query->Y(), query->Z()};

    std::vector<size_t> ret_index(10);
    std::vector<double> out_dist_sqr(10);
    nanoflann::KNNResultSet<double> resultSet(10);
    resultSet.init(ret_index.data(), out_dist_sqr.data());

    nanoflann::SearchParameters params;
    params.sorted = true;
    m_index.findNeighbors(resultSet, queryPoint, params);

    size_t found = resultSet.size();
    for (size_t i = 0; i < found; ++i) {
      Parton* p = m_cloud.get(ret_index[i]);
      if (p->IsUsed()) continue;
      if (p->GetBaryonNumber() * query->GetBaryonNumber() < 0) { return p; }
    }
    return nullptr;
  }

  std::vector<Parton*> FindNearestSame(const Parton* query, size_t k) const {
    double queryPoint[3] = {query->X(), query->Y(), query->Z()};

    std::vector<size_t> ret_index(k * 5);  // 查多点，留空间排除
    std::vector<double> out_dist_sqr(k * 5);
    nanoflann::KNNResultSet<double> resultSet(k * 5);
    resultSet.init(ret_index.data(), out_dist_sqr.data());

    nanoflann::SearchParameters params;
    params.sorted = true;
    m_index.findNeighbors(resultSet, queryPoint, params);

    std::vector<Parton*> result;
    size_t found = resultSet.size();
    for (size_t i = 0; i < found && result.size() < k; ++i) {
      Parton* p = m_cloud.get(ret_index[i]);
      if (p == query) continue;  // ❗️排除自己
      if (p->IsUsed()) continue;
      if (p->GetBaryonNumber() * query->GetBaryonNumber() > 0) { result.push_back(p); }
    }
    return result;
  }

 private:
  PartonCloud m_cloud;
  KDTree m_index;
};

#endif  // PARTON_KDTREE_H
