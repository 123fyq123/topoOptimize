#ifndef _DISTANCE_MAP_
#define _DISTANCE_MAP_

#include "common.h"

struct Point
{
	int dx, dy;

	int DistSq() const { return dx * dx + dy * dy; }
};
struct Cmp_map_br
{
	bool operator () (const Eigen::Vector2i  &a, const Eigen::Vector2i  &b)const
	{
		if (a(0) != b(0))
			return a(0) < b(0);
		else return a(1) < b(1);
	}
};
struct CloestPoint {
	
	//std::set<Eigen::Vector2i, Cmp_map_br> closet_queue;
	std::vector<Eigen::Vector2i> cloestp;
};
struct Grid
{
	std::vector<std::vector<Point>> grid;
	Eigen::MatrixXd dis;
	//Point grid[HEIGHT][WIDTH];
	//std::vector<std::vector<CloestPoint>> cp;
	std::vector<std::vector<CloestPoint>> cp;
};

class DistanceMap
{
public:
	DistanceMap(uint32_t width, uint32_t height,const Eigen::MatrixXi& init_set) :width(width), height(height) {
		grid1.grid.resize(height);
		grid1.cp.resize(height);
		grid1.dis.resize(height,width);
		grid2.grid.resize(height);
		grid2.cp.resize(height);
		grid2.dis.resize(height, width);

		for (uint32_t i = 0; i < height; ++i) {
			grid1.grid[i].resize(width);
			grid1.cp[i].resize(width);
			grid2.cp[i].resize(width);
			grid2.grid[i].resize(width);
		}

		//load to grid1 and grid2

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				
				// Points inside get marked with a dx/dy of zero.
				// Points outside get marked with an infinitely large distance.
				if (init_set(y,x)>=0) {
					Put(grid1, x, y, inside);
					
					Put(grid2, x, y, empty_);
					
				}
				else {
					Put(grid2, x, y, inside);
					
					Put(grid1, x, y, empty_);
					
				}
				//grid1.cp[y][x].closet_queue.push_back(Eigen::Vector2i(y, x));
				//grid2.cp[y][x].closet_queue.push_back(Eigen::Vector2i(y, x));
			}
		}
	};
	~DistanceMap() {};

	inline Point Get(Grid& g, int x, int y)
	{
		// OPTIMIZATION: you can skip the edge check code if you make your grid
		// have a 1-pixel gutter.
		if (x >= 0 && y >= 0 && x < width && y < height)
			return g.grid[y][x];
		else
			return empty_;
	}


	inline void Put(Grid& g, int x, int y, const Point& p)
	{
		g.grid[y][x] = p;
	}


	inline void Compare(Grid& g, Point& p, int x, int y, int offsetx, int offsety)
	{
		Point other = Get(g, x + offsetx, y + offsety);
		other.dx += offsetx;
		other.dy += offsety;
		int disc = 0,disc2 = 0;
		//update distance and cloest points
		disc = other.DistSq();
		disc2 = p.DistSq();
		

		if (disc < disc2) {
			p = other;
			/*
			if ((disc > 10) || (disc == 0)) return;
			//std::cout << y << " " << x << " " << y + offsety << " " << x + offsetx << "\n";
			if (g.cp[y + offsety][x + offsetx].closet_queue.size() == 0) {//init

				std::set<Eigen::Vector2i, Cmp_map_br>().swap(g.cp[y][x].closet_queue);
				g.cp[y][x].closet_queue.insert(Eigen::Vector2i(y + offsety, x + offsetx));
					
			}
			else {//reset
				std::set<Eigen::Vector2i, Cmp_map_br>().swap(g.cp[y][x].closet_queue);
				for (auto i = g.cp[y + offsety][x + offsetx].closet_queue.begin(); i != g.cp[y + offsety][x + offsetx].closet_queue.end(); ++i) {
					g.cp[y][x].closet_queue.insert(*i);
				}
			}
			*/
			
			
		}
		else if ( disc == disc2) {
			//距离相等，记录多个
			/*
			if ((disc > 10) || (disc == 0)) return;
			//std::cout << y << " " << x << " " << y + offsety << " " << x + offsetx << "\n";
			for (auto i = g.cp[y + offsety][x + offsetx].closet_queue.begin(); i != g.cp[y + offsety][x + offsetx].closet_queue.end(); ++i) {
				g.cp[y][x].closet_queue.insert(*i);
			}
			*/
		}
	}



	inline void GenerateSDF(Grid& g)
	{
		// Pass 0
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				Point p = Get(g, x, y);
				Compare(g, p, x, y, -1, 0);
				Compare(g, p, x, y, 0, -1);
				Compare(g, p, x, y, -1, -1);
				Compare(g, p, x, y, 1, -1);
				Put(g, x, y, p);
			}

			for (int x = width - 1; x >= 0; x--)
			{
				Point p = Get(g, x, y);
				Compare(g, p, x, y, 1, 0);
				Put(g, x, y, p);
			}
		}

		// Pass 1
		for (int y = height - 1; y >= 0; y--)
		{
			for (int x = width - 1; x >= 0; x--)
			{
				Point p = Get(g, x, y);
				Compare(g, p, x, y, 1, 0);
				Compare(g, p, x, y, 0, 1);
				Compare(g, p, x, y, -1, 1);
				Compare(g, p, x, y, 1, 1);
				Put(g, x, y, p);
			}

			for (int x = 0; x < width; x++)
			{
				Point p = Get(g, x, y);
				Compare(g, p, x, y, -1, 0);
				Put(g, x, y, p);
			}
		}

		//collect cloest point
		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				//std::cout << "y:" << y << " x:" << x << "\n";
				//if((g.grid[y][x].dy!=0)&&(g.grid[y][x].dx!=0))
				g.cp[y][x].cloestp.push_back(Eigen::Vector2i(y+g.grid[y][x].dy, x + g.grid[y][x].dx));
			}
		}
	}

	inline void Generate_Bdy_inner(bool inner, bool bdy) {
		if (inner) {
			GenerateSDF(grid1);
		}
		if (bdy) {
			GenerateSDF(grid2);
		}
		grid_value.resize(height,width);
		double dist1, dist2, dist;
		for (uint32_t row = 0; row < height; ++row) {
			for (uint32_t col = 0; col < width; ++col) {
				dist1 = (std::sqrt((double)Get(grid1, col, row).DistSq()));
				dist2 = (std::sqrt((double)Get(grid2, col, row).DistSq()));
				dist = dist1 - dist2;
				grid1.dis(row, col) = dist;
				
				grid2.dis(row, col) = -dist;
			}
		}
		

	};

	Grid grid1, grid2;
	Eigen::MatrixXd grid_value;
private:
	uint32_t width;
	uint32_t height;
	
	const Point inside = { 0,0 };
	const Point empty_ = { 9999,9999 };
	const int maxdis = 9999 * 9999 + 9999 * 9999;
};




#endif // !_DISTANCE_MAP_

