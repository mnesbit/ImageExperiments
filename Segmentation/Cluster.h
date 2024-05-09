#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include <vector>
#include <memory>
#include "../ImageHelper/inc/image.h"

namespace cluster {

	struct ClusterItem;
	struct Edge;

	class Clusters
	{
	public:
		Clusters(const img::image<img::rgb>* input);
		~Clusters();
		void Process(int rounds);
		std::unique_ptr<img::image<img::rgb> > Averages() const;
		std::unique_ptr<img::image<img::rgb> > Borders() const;
	private:
		size_t m_height;
		size_t m_width;
		size_t m_clusterCount;
		size_t m_level;
		size_t m_round;
		std::vector<Edge> m_edges;
		std::vector<ClusterItem> m_clusters;

		void InitClusters(const img::image<img::rgb>* input);
		void InitEdges();
		void CalcGradient();
		size_t Find(size_t a) const;
		size_t Next(size_t a) const;
		void Mark(size_t a, int value);
		int GetMark(size_t a) const;
		void ClearMarks();
		void Join(const Edge& edge);
		void SortEdges();
	};

}

#endif //_CLUSTER_H_