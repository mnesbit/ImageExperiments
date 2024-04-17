#include "cluster.h"
#include "../ImageHelper/inc/misc.h"
#include <algorithm>
#include <iostream>
#include <format>
#include <unordered_map>

namespace cluster {
	//Farid, H., & Simoncelli, E. P. (2004). Differentiation of discrete multidimensional signals. IEEE Transactions on image processing, 13(4), 496-508.
	static const double g_diff_3[3] { -0.425287,0.000000,0.425287 }; // length 3
	static const double g_interp_3[3] { 0.229879,0.540242,0.229879 };
	static const double g_diff_5[5] { -0.109604,-0.276691,0.000000,0.276691,0.109604 }; // length 5
	static const double g_interp_5[5] { 0.037659,0.249153,0.426375,0.249153,0.037659 };
	static const double g_diff_7[7] { -0.015964,-0.121482,-0.193357,0.000000,0.193357,0.121482,0.015964 }; // length 7
	static const double g_interp_7[7] { 0.003992,0.067088,0.246217,0.365406,0.246217,0.067088,0.003992 };
	static const double g_diff_9[9] { -0.003059,-0.035187,-0.118739,-0.143928,0.000000,0.143928,0.118739,0.035187,0.003059 }; // length 9
	static const double g_interp_9[9] { 0.000721,0.015486,0.090341,0.234494,0.317916,0.234494,0.090341,0.015486,0.000721 };

	struct ClusterItem
	{
		size_t next;
		size_t prev;
		size_t rank;
		mutable size_t parent;
		int mark;
		img::rgb orig;
		double origy;
		double sumr;
		double sumg;
		double sumb;
		double gradx;
		double grady;
		size_t size;
	};

	struct Edge
	{
		size_t a;
		size_t b;
		int dir;
		double grad;
		double weight;
	};

	inline bool operator<(const Edge& a, const Edge& b)
	{
		return (a.weight < b.weight);
	}

	void Clusters::InitEdges()
	{
		m_edges.reserve(4 * m_height * m_width);
		for (size_t y = 0; y < m_height; y++)
		{
			for (size_t x = 0; x < m_width; x++)
			{
				const ClusterItem& clusta = m_clusters[y * m_width + x];
				if (x < m_width - 1)
				{
					const ClusterItem& clustb = m_clusters[y * m_width + (x + 1)];
					double totx = square(1.0 * (clusta.gradx + clustb.gradx)); // summed gradient in x direction. The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					double toty = square(1.0 * (clusta.grady + clustb.grady)); // summed gradient in y direction The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					m_edges.push_back(
						Edge{
						.a = y * m_width + x ,
						.b = y * m_width + (x + 1) ,
						.dir = 0,
						.grad = sqrt(totx + toty)
						}
					);
				}
				if (y < m_height - 1)
				{
					const ClusterItem& clustb = m_clusters[(y + 1) * m_width + x];
					double totx = square(1.0 * (clusta.gradx + clustb.gradx)); // summed gradient in x direction. The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					double toty = square(1.0 * (clusta.grady + clustb.grady)); // summed gradient in y direction The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					m_edges.push_back(
						Edge{
						.a = y * m_width + x ,
						.b = (y + 1) * m_width + x ,
						.dir = 1,
						.grad = sqrt(totx + toty)
						}
					);
				}
				if ((x < m_width - 1) && (y < m_height - 1))
				{
					const ClusterItem& clustb = m_clusters[(y + 1) * m_width + (x + 1)];
					double totx = square(1.0 * (clusta.gradx + clustb.gradx)); // summed gradient in x direction. The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					double toty = square(1.0 * (clusta.grady + clustb.grady)); // summed gradient in y direction The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					m_edges.push_back(
						Edge{
						.a = y * m_width + x,
						.b = (y + 1) * m_width + (x + 1),
						.dir = 2,
						.grad = sqrt(totx + toty)
						}
					);
				}
				if ((x < m_width - 1) && (y > 0))
				{
					const ClusterItem& clustb = m_clusters[(y - 1) * m_width + (x + 1)];
					double totx = square(1.0 * (clusta.gradx + clustb.gradx)); // summed gradient in x direction. The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					double toty = square(1.0 * (clusta.grady + clustb.grady)); // summed gradient in y direction The weighting factor of 1.0 rather than 0.5 seems a good balance between edge and patch energy
					m_edges.push_back(
						Edge{
						.a = y * m_width + x ,
						.b = (y - 1) * m_width + (x + 1) ,
						.dir = 3,
						.grad = sqrt(totx + toty)
						}
					);
				}
			}
		}
	}

	void Clusters::CalcGradient()
	{
		for (size_t x = 0; x < m_width; x++)
		{
			for (size_t y = 0; y < m_height; y++)
			{
				double totxr = 0.0;
				double totyr = 0.0;
				double totxg = 0.0;
				double totyg = 0.0;
				double totxb = 0.0;
				double totyb = 0.0;
				for (int dx = -4; dx <= +4; dx++)
				{
					int u = std::clamp(static_cast<int>(x) + dx, 0, static_cast<int>(m_width - 1));
					for (int dy = -4; dy <= +4; dy++)
					{
						int v = std::clamp(static_cast<int>(y) + dy, 0, static_cast<int>(m_height - 1));
						img::rgb pt = m_clusters[v * m_width + u].orig;
						img::yuv pt2 = img::YUVFromRGB(pt);
						if (dx >= -2
							&& dx <= 2
							&& dy >= -2
							&& dy <= 2) // use 5 x 5 differentiator for Y channel
						{
							double coeff = g_diff_5[dx + 2] * g_interp_5[dy + 2];
							double coeff2 = g_diff_5[dy + 2] * g_interp_5[dx + 2];
							totxr += coeff * pt2.y;
							totyr += coeff2 * pt2.y;
						}
						if (dx >= -3
							&& dx <= 3
							&& dy >= -3
							&& dy <= 3) // use 7 x 7 differentiator for U channel
						{
							double coeff = g_diff_7[dx + 3] * g_interp_7[dy + 3];
							double coeff2 = g_diff_7[dy + 3] * g_interp_7[dx + 3];
							totxg += coeff * pt2.u;
							totyg += coeff2 * pt2.u;
						}
						{
							// use 9 x 9 differentiator for V channel
							double coeff = g_diff_9[dx + 4] * g_interp_9[dy + 4];
							double coeff2 = g_diff_9[dy + 4] * g_interp_9[dx + 4];
							totxb += coeff * pt2.v;
							totyb += coeff2 * pt2.v;
						}
					}
				}
				m_clusters[y * m_width + x].gradx = abs(totxr) + abs(totxg) + abs(totxb);
				m_clusters[y * m_width + x].grady = abs(totyr) + abs(totyg) + abs(totyb);
			}
		}
	}

	void Clusters::InitClusters(const img::image<img::rgb>* input) {
		m_clusterCount = m_width * m_height;
		m_clusters.resize(m_clusterCount);
		for (size_t y = 0; y < m_height; y++)
		{
			for (size_t x = 0; x < m_width; x++)
			{
				img::rgb point = imRef(input, x, y);
				m_clusters[y * m_width + x].orig = point;
				double pt = YFromRGB(point);
				m_clusters[y * m_width + x].origy = pt;
			}
		}
		for (size_t i = 0; i < m_clusterCount; i++)
		{
			m_clusters[i].rank = 0;
			m_clusters[i].mark = -1;
			m_clusters[i].parent = i; // initially every pixel is a cluster
			m_clusters[i].prev = i - 1;
			m_clusters[i].next = i + 1;
			m_clusters[i].size = 1;
			m_clusters[i].sumr = m_clusters[i].orig.r;
			m_clusters[i].sumg = m_clusters[i].orig.g;
			m_clusters[i].sumb = m_clusters[i].orig.b;
		}
		m_clusters[0].prev = m_clusterCount - 1;
		m_clusters[m_clusterCount - 1].next = 0;
	}

	Clusters::Clusters(const img::image<img::rgb>* input)
		: m_height(input->height())
		, m_width(input->width())
		, m_level(0)
		, m_clusterCount(0)
		, m_round(0)
	{
		InitClusters(input);
		CalcGradient();
		InitEdges();
	}

	Clusters::~Clusters() {
	}

	size_t Clusters::Find(size_t a) const
	{
		size_t cluster = a;
		while (cluster != m_clusters[cluster].parent) // find cluster root
		{
			cluster = m_clusters[cluster].parent;
		}
		m_clusters[a].parent = cluster; // prefix compression
		return cluster;
	}

	size_t Clusters::Next(size_t a) const
	{
		return m_clusters[a].next;
	}

	void Clusters::Mark(size_t a, int value)
	{
		m_clusters[Find(a)].mark = value;
	}

	int Clusters::GetMark(size_t a) const
	{
		return m_clusters[Find(a)].mark;
	}

	void Clusters::ClearMarks()
	{
		size_t first = Find(0);
		Mark(first, -1);
		for (size_t curr = Next(first); curr != first; curr = Next(curr))
		{
			Mark(curr, -1);
		}
	}

	void Clusters::Join(const Edge& edge)
	{
		size_t a = Find(edge.a);
		size_t b = Find(edge.b);
		if (m_clusters[a].rank > m_clusters[b].rank)
		{
			m_clusters[b].parent = a;
			m_clusters[m_clusters[b].prev].next = m_clusters[b].next;
			m_clusters[m_clusters[b].next].prev = m_clusters[b].prev;
			m_clusters[a].size += m_clusters[b].size;
			m_clusters[a].sumr += m_clusters[b].sumr;
			m_clusters[a].sumg += m_clusters[b].sumg;
			m_clusters[a].sumb += m_clusters[b].sumb;
		}
		else
		{
			m_clusters[a].parent = b;
			m_clusters[m_clusters[a].prev].next = m_clusters[a].next;
			m_clusters[m_clusters[a].next].prev = m_clusters[a].prev;
			m_clusters[b].size += m_clusters[a].size;
			m_clusters[b].sumr += m_clusters[a].sumr;
			m_clusters[b].sumg += m_clusters[a].sumg;
			m_clusters[b].sumb += m_clusters[a].sumb;
			if (m_clusters[a].rank == m_clusters[b].rank)
			{
				m_clusters[b].rank++;
				if (m_clusters[b].rank > m_level)
				{
					m_level = m_clusters[b].rank;
				}
			}
		}
		m_clusterCount--;
	}

	void Clusters::SortEdges() {
		std::erase_if(m_edges, [&](const Edge& edge) {
			size_t a = Find(edge.a);
			size_t b = Find(edge.b);
			return a == b;
		});
		for (auto& edge : m_edges) {
			size_t a = Find(edge.a);
			size_t b = Find(edge.b);
			double sza = (double)m_clusters[a].size;
			double szb = (double)m_clusters[b].size;
			img::yuv acol = img::YUVFromRGB((m_clusters[a].sumr / sza), (m_clusters[a].sumg / sza), (m_clusters[a].sumb / sza));
			img::yuv bcol = img::YUVFromRGB((m_clusters[b].sumr / szb), (m_clusters[b].sumg / szb), (m_clusters[b].sumb / szb));
			double tot = square(acol.y - bcol.y) + square(acol.u - bcol.u) + square(acol.v - bcol.v); //squared colour difference of adjacent regions
			edge.weight = edge.grad + sqrt(tot); // by summing gradient and region difference the gradient helps form sensible groups in early rounds and layter rounds use averaged colour
		}
		std::sort(m_edges.begin(), m_edges.end()); // sort by weight so that strongest differences are merged first
	}

	void Clusters::Process(int rounds)
	{
		SortEdges();
		// segment
		int joins;
		double n = static_cast<double>(m_clusterCount);
		std::cout << std::format("Initial components {}\n", m_clusterCount) << std::endl;
		for (int rep = 0;rep < rounds;++rep)
		{
			ClearMarks();
			joins = 0;
			if (m_edges.empty()) break;
			for (auto& edge: m_edges)
			{
				size_t a = Find(edge.a);
				size_t b = Find(edge.b);
				if (a != b)
				{
					if ((GetMark(a) == -1) || (GetMark(b) == -1)) // we require at least one of the mergees to be untouched in this round.
					{
						Join(edge);
						joins++;
						Mark(a, 1);
					}
				}
			}
			SortEdges();
			std::cout << std::format("edges {} components {} ratio {} joins {}", m_edges.size(), m_clusterCount, (double)m_clusterCount / n, joins) << std::endl;
			n = static_cast<double>(m_clusterCount);
		}
	}

	img::image<img::rgb>* Clusters::Averages() const
	{
		img::image<img::rgb>* output = new img::image<img::rgb>(m_width, m_height);
		for (size_t y = 0; y < m_height; y++)
		{
			for (size_t x = 0; x < m_width; x++)
			{
				size_t comp = Find(y * m_width + x);
				double size = static_cast<double>(m_clusters[comp].size);
				imRef(output, x, y).r = static_cast<img::uchar>(m_clusters[comp].sumr / size);
				imRef(output, x, y).g = static_cast<img::uchar>(m_clusters[comp].sumg / size);
				imRef(output, x, y).b = static_cast<img::uchar>(m_clusters[comp].sumb / size);
			}
		}
		return output;
	}

	img::image<img::rgb>* Clusters::Borders() const
	{
		img::image<img::rgb>* output = new img::image<img::rgb>(m_width, m_height, true);
		std::unordered_map<size_t, img::rgb>  colors;
		for (auto& edge : m_edges)
		{
			size_t a = Find(edge.a);
			size_t b = Find(edge.b);
			if (!colors.contains(a)) {
				img::rgb rancolour{
					.r = (img::uchar)rand(),
					.g = (img::uchar)rand(),
					.b = (img::uchar)rand()
				};
				colors.insert(std::make_pair(a, rancolour));
			}
			if (!colors.contains(b)) {
				img::rgb rancolour{
					.r = (img::uchar)rand(),
					.g = (img::uchar)rand(),
					.b = (img::uchar)rand()
				};
				colors.insert(std::make_pair(b, rancolour));
			}
			if (a != b)
			{
				size_t ay = edge.a / m_width;
				size_t ax = edge.a - (m_width * ay);
				img::rgb cola = colors[a];
				imRef(output, ax, ay).r = cola.r;
				imRef(output, ax, ay).g = cola.g;
				imRef(output, ax, ay).b = cola.b;
				size_t by = edge.b / m_width;
				size_t bx = edge.b - (m_width * by);
				img::rgb colb = colors[b];
				imRef(output, bx, by).r = colb.r;
				imRef(output, bx, by).g = colb.g;
				imRef(output, bx, by).b = colb.b;
			}
		}
		return output;
	}
}