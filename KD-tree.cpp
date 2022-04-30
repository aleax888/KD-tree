#include <iostream>
#include <utility>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <climits>
#include <random>
using namespace std;

/*generates a random number between two numbers that it receives as a 
parameter, an upper limit and a lower limit*/
int generate_ramdon_number(int dawn, int up)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<> distrib(dawn, up);
	return distrib(gen);
}

struct point
{
	//attributes
	int dimension;
	double *coordinates;

	//constructor
	point() // void constructor
	{
		int dimension = 0;
		coordinates = new double[0];
	}

	point(int d) //semi random constructor
	{
		dimension = d;
		coordinates = new double[d];
		generate_ramdon_array(coordinates, dimension);
	}

	point(int d, double *coord) // constructor
	{
		dimension = d;
		coordinates = new double[dimension];
		for (int i = 0; i < dimension; ++i)
			coordinates[i] = coord[i];
	}

	point(const point &p) // copy constructor
	{
		dimension = p.dimension;
		coordinates = new double[dimension];
		for (int i = 0; i < dimension; ++i)
			coordinates[i] = p.coordinates[i];
	}

	//destructor
	/*~point()
	{
		delete[] coordinates;
	}*/

	//methods
	void copy(point p) // copies the coordinates of a point received as a parameter
	{
		dimension = p.dimension;
		coordinates = new double[dimension];
		for (int i = 0; i < dimension; ++i)
			coordinates[i] = p.coordinates[i];
	}

	bool equal(point p) // return true if p1 and p2 are the same, else return false
	{
		if (this->dimension != p.dimension)
			return false;
		for (int i = 0; i < this->dimension; ++i)
			if (this->coordinates[i] != p.coordinates[i])
				return false;
		return true;
	}

	bool distinct(point p) // return false if p1 and p2 are the same, else return true
	{
		return !equal(p);
	}

	double euclidean_distance(point p) // return distance between two points
	{
		double ans = 0;
		for (int i = 0; i < p.dimension; ++i)
		{
			ans += pow(this->coordinates[i] - p.coordinates[i], 2);
		}
		return sqrt(ans);
	}

	void printPoint() // print coordinates of the point
	{
		cout << "Point: ";
		for (int i = 0; i < dimension; ++i)
			cout << coordinates[i] << " ";
		cout << endl;
	}

	void generate_ramdon_array(double* arr, int size) // it is used for the constructor that does not receive an array, it generates an array of random coordinates
	{
		random_device rd;
		mt19937 gen(rd());
		uniform_int_distribution<> distrib(1, 1000);
		for (int i = 0; i < size; ++i)
			arr[i] = distrib(gen);
	}
};

struct kd_node
{
	//attributes
	point datum;
	kd_node *child[2] = {};

	//constructor
	kd_node(point p)
	{
		datum.copy(p);
		child[0] = nullptr;
		child[1] = nullptr;
	}
};

// function to change the sort order of the sort algorithm
bool sortbysec(const pair<kd_node*, double>& a, const pair<kd_node*, double>& b)
{
	return (a.second < b.second);
}

struct kd_tree
{
	//attributes
	kd_node *root;
	int dimension;

	//constructor
	kd_tree(int d)
	{
		root = nullptr;
		dimension = d;
	}

	//methods

	//visualization
	void in_orden(kd_node *kdn)
	{
		if (kdn)
		{
			in_orden(kdn->child[0]);
			kdn->datum.printPoint();
			in_orden(kdn->child[1]);
		}
	}

	void pre_orden(kd_node* kdn)
	{
		if (kdn)
		{
			kdn->datum.printPoint();
			in_orden(kdn->child[0]);
			in_orden(kdn->child[1]);
		}
	}

	void post_orden(kd_node* kdn)
	{
		if (kdn)
		{
			in_orden(kdn->child[0]);
			in_orden(kdn->child[1]);
			kdn->datum.printPoint();
		}
	}

	
	bool serch(point p, kd_node **&finder, int coord = 0) // find a pointer by a point
	{
		coord = coord % dimension;

		if (*finder == nullptr)
			return false;
		
		if ((*finder)->datum.equal(p))
			return true;
		
		if (p.coordinates[coord] >= (*finder)->datum.coordinates[coord])
		{
			finder = &(*finder)->child[1];
			serch(p, finder, coord + 1);
		}
		else
		{
			finder = &(*finder)->child[0];
			serch(p, finder, coord + 1);
		}
	}

	void insert(point p) // first use function serch to find pointer to insert a new point
	{
		kd_node** finder = &root;
		if (!serch(p, finder))
			*finder = new kd_node(p);
	}

	void erase(point p) // erase like binary tree
	{
		kd_node** finder = &root;
		if (serch(p, finder))
		{
			kd_node** aux = finder;
			if ((*aux)->child[0])
			{
				aux = &(*aux)->child[0];
				
				while ((*aux)->child[0] or (*aux)->child[1])
				{
					while ((*aux)->child[1])
						aux = &(*aux)->child[1];
					if ((*aux)->child[0])
						aux = &(*aux)->child[0];
				}
				
				(*finder)->datum.copy((*aux)->datum);
				*aux = nullptr;
			}

			else if ((*aux)->child[1])
			{
				aux = &(*aux)->child[1];
				while ((*aux)->child[0] or (*aux)->child[1])
				{
					while ((*aux)->child[0])
						aux = &(*aux)->child[0];
					if ((*aux)->child[1])
						aux = &(*aux)->child[1];
				}
				(*finder)->datum.copy((*aux)->datum);
				*aux = nullptr;
			}

			else
				*finder = nullptr;
		}
	}

	void nearest_neighbor(kd_node *target, kd_node *current_node, kd_node *&nearest_neighbors_candidate, int depth = 0,	double best_distance = LLONG_MAX)
	{
		if (!current_node)
			return;
		
		if (current_node->datum.euclidean_distance(target->datum) < best_distance)
		{
			best_distance = current_node->datum.euclidean_distance(target->datum);
			nearest_neighbors_candidate = current_node;
		}

		int axis = depth % dimension;
		bool right = false;

		if (target->datum.coordinates[axis] > current_node->datum.coordinates[axis]) 
		{	
			right = true;
			nearest_neighbor(target, current_node->child[1], nearest_neighbors_candidate, depth+1, best_distance);
		}
		
		else 
		{
			right = false;
			nearest_neighbor(target, current_node->child[0], nearest_neighbors_candidate, depth + 1, best_distance);
		}

		if (fabs(current_node->datum.coordinates[axis] - target->datum.coordinates[axis]) < best_distance) 
		{
			if (right) 
				nearest_neighbor(target, current_node->child[0], nearest_neighbors_candidate, depth + 1, best_distance);
			else 
				nearest_neighbor(target, current_node->child[1], nearest_neighbors_candidate, depth + 1, best_distance);
		}
	}
	
	bool compare(vector<pair<kd_node*, double>>& k_nearest_neighbors_candidates, double d)
	{
		for (int i = 0, s = k_nearest_neighbors_candidates.size(); i < s; ++i)
		{
			if (d < k_nearest_neighbors_candidates[i].second)
				return true;
		}
		return false;
	}
	void k_nearest_neighbor(kd_node* target, int k, kd_node* current_node,	vector<pair<kd_node*,double>> &k_nearest_neighbors_candidates, int depth = 0)
	{
		if (!current_node)
			return;

		k_nearest_neighbors_candidates.push_back(pair<kd_node*, double>{current_node, current_node->datum.euclidean_distance(target->datum)});
		sort(k_nearest_neighbors_candidates.begin(), k_nearest_neighbors_candidates.end(), sortbysec);
		if (k_nearest_neighbors_candidates.size() > k) k_nearest_neighbors_candidates.pop_back();
		
		int axis = depth % dimension;
		bool right = false;
		if (target->datum.coordinates[axis] < current_node->datum.coordinates[axis]) 
		{
			right = true;
			k_nearest_neighbor(target, k, current_node->child[1], k_nearest_neighbors_candidates, depth + 1);
		}
		
		else 
		{
			right = false;
			k_nearest_neighbor(target, k, current_node->child[0], k_nearest_neighbors_candidates, depth + 1);
		}

		if (compare(k_nearest_neighbors_candidates, fabs(current_node->datum.coordinates[axis] - target->datum.coordinates[axis])))
		{
			if (right) {
				k_nearest_neighbor(target, k, current_node->child[0], k_nearest_neighbors_candidates, depth + 1);
			}
			else {
				k_nearest_neighbor(target, k, current_node->child[1], k_nearest_neighbors_candidates, depth + 1);
			}
		}
	}

	void range_query(kd_node* target, double query_distance, kd_node* current_node,	vector<pair<kd_node*, double>>& k_nearest_neighbors_candidates, int depth = 0)
	{
		if (!current_node)
			return;

		if (current_node->datum.euclidean_distance(target->datum) <= query_distance)
		{
			k_nearest_neighbors_candidates.push_back(pair<kd_node*, double>{current_node, current_node->datum.euclidean_distance(target->datum)});
			sort(k_nearest_neighbors_candidates.begin(), k_nearest_neighbors_candidates.end(), sortbysec);
		}
		
		int axis = depth % dimension;
		bool right = false;
		if (target->datum.coordinates[axis] > current_node->datum.coordinates[axis])
		{
			right = true;
			range_query(target, query_distance, current_node->child[1], k_nearest_neighbors_candidates, depth + 1);
		}

		else
		{
			right = false;
			range_query(target, query_distance, current_node->child[0], k_nearest_neighbors_candidates, depth + 1);
		}

		if (fabs(current_node->datum.coordinates[axis] - target->datum.coordinates[axis]) < query_distance)
		{
			if (right) {
				range_query(target, query_distance, current_node->child[0], k_nearest_neighbors_candidates, depth + 1);
			}
			else {
				range_query(target, query_distance, current_node->child[1], k_nearest_neighbors_candidates, depth + 1);
			}
		}
	}

}; 

//test by a example
void test_insert()
{
	double arr1[2] = { 40, 45 };
	point p1(2, arr1);
	double arr2[2] = { 15, 70 };
	point p2(2, arr2);
	double arr3[2] = { 70, 10 };
	point p3(2, arr3);
	double arr4[2] = { 69, 50 };
	point p4(2, arr4);
	double arr5[2] = { 66, 85 };
	point p5(2, arr5);
	double arr6[2] = { 85, 90 };
	point p6(2, arr6);

	kd_tree tree(2);
	tree.insert(p1);
	tree.insert(p2);
	tree.insert(p3);
	tree.insert(p4);
	tree.insert(p5);
	tree.insert(p6);

	tree.in_orden(tree.root);
}
void test_erase()
{
	double arr1[2] = { 40, 45 };
	point p1(2, arr1);
	double arr2[2] = { 15, 70 };
	point p2(2, arr2);
	double arr3[2] = { 70, 10 };
	point p3(2, arr3);
	double arr4[2] = { 69, 50 };
	point p4(2, arr4);
	double arr5[2] = { 66, 85 };
	point p5(2, arr5);
	double arr6[2] = { 85, 90 };
	point p6(2, arr6);

	kd_tree tree(2);
	tree.insert(p1);
	tree.insert(p2);
	tree.insert(p3);
	tree.insert(p4);
	tree.insert(p5);
	tree.insert(p6);

	tree.in_orden(tree.root);
	tree.erase(p1);
	cout << "erase ---" << endl;
	tree.in_orden(tree.root);
	tree.erase(p2);
	cout << "erase ---" << endl;
	tree.in_orden(tree.root);
}
void test_serch()
{
	double arr1[2] = { 40, 45 };
	point p1(2, arr1);
	double arr2[2] = { 15, 70 };
	point p2(2, arr2);
	double arr3[2] = { 70, 10 };
	point p3(2, arr3);
	double arr4[2] = { 69, 50 };
	point p4(2, arr4);
	double arr5[2] = { 66, 85 };
	point p5(2, arr5);
	double arr6[2] = { 85, 90 };
	point p6(2, arr6);


	kd_tree tree(2);
	tree.insert(p1);
	tree.insert(p2);
	tree.insert(p3);
	tree.insert(p4);
	tree.insert(p5);
	tree.insert(p6);

	kd_node **finder = &tree.root;
	double arr7[2] = { 0, 0 };
	point p7(2, arr7);
	cout << "serch: " << tree.serch(p6, finder) << ", ";
	cout << (*finder)->datum.coordinates[0] << " " << (*finder)->datum.coordinates[1];
}
void test_nn()
{
	double arr1[2] = { 40, 45 };
	point p1(2, arr1);
	double arr2[2] = { 15, 70 };
	point p2(2, arr2);
	double arr3[2] = { 70, 10 };
	point p3(2, arr3);
	double arr4[2] = { 69, 50 };
	point p4(2, arr4);
	double arr5[2] = { 66, 85 };
	point p5(2, arr5);
	double arr6[2] = { 85, 90 };
	point p6(2, arr6);

	double arr7[2] = { 80, 80 };
	point p7(2, arr7);
	kd_node *target = new kd_node(p7), *nn = nullptr;

	kd_tree tree(2);
	tree.insert(p1);
	tree.insert(p2);
	tree.insert(p3);
	tree.insert(p4);
	tree.insert(p5);
	tree.insert(p6);

	//tree.in_orden(tree.root);
	tree.nearest_neighbor(target, tree.root, nn);
	cout << "target           -> "; target->datum.printPoint();
	cout << "nearest neighbor -> "; nn->datum.printPoint();
}
void test_knn()
{
	double arr1[2] = { 40, 45 };
	point p1(2, arr1);
	double arr2[2] = { 15, 70 };
	point p2(2, arr2);
	double arr3[2] = { 70, 10 };
	point p3(2, arr3);
	double arr4[2] = { 69, 50 };
	point p4(2, arr4);
	double arr5[2] = { 66, 85 };
	point p5(2, arr5);
	double arr6[2] = { 85, 90 };
	point p6(2, arr6);

	double arr7[2] = { 80, 80 };
	point p7(2, arr7);
	kd_node* target = new kd_node(p7);
	vector<pair<kd_node*, double>> knn;

	kd_tree tree(2);
	tree.insert(p1);
	tree.insert(p2);
	tree.insert(p3);
	tree.insert(p4);
	tree.insert(p5);
	tree.insert(p6);

	//tree.in_orden(tree.root);
	tree.k_nearest_neighbor(target, 3, tree.root, knn);
	cout << "target             -> "; target->datum.printPoint();
	cout << "k nearest neighbor: " << endl;
	for (int i = 0; i < knn.size(); ++i)
	{
		knn[i].first->datum.printPoint(); cout << "-> distance: " << knn[i].second << endl;
	}
}
void test_rq()
{
	double arr1[2] = { 40, 45 };
	point p1(2, arr1);
	double arr2[2] = { 15, 70 };
	point p2(2, arr2);
	double arr3[2] = { 70, 10 };
	point p3(2, arr3);
	double arr4[2] = { 69, 50 };
	point p4(2, arr4);
	double arr5[2] = { 66, 85 };
	point p5(2, arr5);
	double arr6[2] = { 85, 90 };
	point p6(2, arr6);

	double arr7[2] = { 80, 80 };
	point p7(2, arr7);
	kd_node* target = new kd_node(p7);
	vector<pair<kd_node*, double>> rq;

	kd_tree tree(2);
	tree.insert(p1);
	tree.insert(p2);
	tree.insert(p3);
	tree.insert(p4);
	tree.insert(p5);
	tree.insert(p6);

	//tree.inorden(tree.root);
	tree.range_query(target, 20, tree.root, rq);
	cout << "target            -> "; target->datum.printPoint();
	cout << "nodes into query radius: " << endl;
	for (int i = 0; i < rq.size(); ++i)
	{
		rq[i].first->datum.printPoint(); cout << "-> distance: " << rq[i].second << endl;
	}
}

// random testing
void test_insert(int size, int d)
{
	kd_tree tree(d);
	vector<point> vec;

	cout << "Number of nodes: " << size << ", Dimension: " << d << endl << endl;
	cout << "Random nodes to insert:" << endl;
	for (int i = 0; i < size; ++i)
	{
		vec.push_back(point(d));
		vec.back().printPoint();
		tree.insert(vec.back());
	}
	cout << endl;

	cout << "In Orden:" << endl; tree.in_orden(tree.root); cout << endl;
	cout << "Pre Orden:" << endl; tree.pre_orden(tree.root); cout << endl;
	cout << "Post Orden:" << endl; tree.post_orden(tree.root); cout << endl;
}

void test_erase(int size, int d, int n_delete)
{
	if (n_delete > size)
		cout << "Number of nodes to erase can't be larger than size" << endl;
	
	else
	{
		kd_tree tree(d);
		vector<point> vec;

		cout << "Number of nodes: " << size << ", Dimension: " << d << endl << endl;
		cout << "Random nodes to insert:" << endl;
		for (int i = 0; i < size; ++i)
		{
			vec.push_back(point(d));
			vec.back().printPoint();
			tree.insert(vec.back());
		}
		cout << endl;

		cout << "In Orden after:" << endl; tree.in_orden(tree.root); cout << endl;

		for (int i = 0, s = vec.size(), random = generate_ramdon_number(0, s - 1); i < n_delete; i++, s--, random = generate_ramdon_number(0, s - 1))
		{
			cout << "Random node to erase: "; vec[random].printPoint(); cout << endl;
			tree.erase(vec[random]);
			vec.erase(vec.begin() + random);
		}

		cout << "In Orden before:" << endl; tree.in_orden(tree.root); cout << endl;
	}
}

void test_serch(int size, int d, bool random)
{
	kd_tree tree(d);
	vector<point> vec;

	cout << "Number of nodes: " << size << ", Dimension: " << d << endl << endl;
	cout << "Random nodes to insert:" << endl;
	for (int i = 0; i < size; ++i)
	{
		vec.push_back(point(d));
		vec.back().printPoint();
		tree.insert(vec.back());
	}
	cout << endl;

	cout << "In Orden:" << endl; tree.in_orden(tree.root); cout << endl << endl;

	if (!random)
	{
		kd_node** finder = &tree.root;
		cout << "Serch: "; vec[0].printPoint(); cout << "Result: " << tree.serch(vec[0], finder) << endl;
		cout << "found: ";  (*finder)->datum.printPoint();
	}

	else
	{
		point p(d);
		kd_node** finder = &tree.root;
		bool result = tree.serch(p, finder);
		cout << "Serch: "; p.printPoint(); cout << "Result: " << result << endl;
		  
		if (tree.serch(vec[0], finder))
		{
			cout << "found: ";
			(*finder)->datum.printPoint();
		}

		else
			cout << "not found" << endl;
	}
}

void test_nn(int size, int d)
{
	kd_tree tree(d);
	vector<point> vec;

	cout << "Number of nodes: " << size << ", Dimension: " << d << endl << endl;
	cout << "Random nodes to insert:" << endl;
	for (int i = 0; i < size; ++i)
	{
		vec.push_back(point(d));
		vec.back().printPoint();
		tree.insert(vec.back());
	}
	cout << endl;

	cout << "In Orden:" << endl; tree.in_orden(tree.root); cout << endl << endl;

	kd_node *target = new kd_node(point(d)), *nn = nullptr;
	tree.nearest_neighbor(target, tree.root, nn);
	cout << "target:           "; target->datum.printPoint();
	cout << "nearest neighbor: "; nn->datum.printPoint();

	cout << endl << "Distances from all nodes to the target: " << endl;

	for (int i = 0; i < size; ++i)
	{
		vec[i].printPoint();
		cout << vec[i].euclidean_distance(target->datum) << endl;
	}
}

void test_knn(int size, int d, int k)
{
	kd_tree tree(d);
	vector<point> vec;

	cout << "Number of nodes: " << size << ", Dimension: " << d << endl << endl;
	cout << "Random nodes to insert:" << endl;
	for (int i = 0; i < size; ++i)
	{
		vec.push_back(point(d));
		vec.back().printPoint();
		tree.insert(vec.back());
	}
	cout << endl;

	cout << "In Orden:" << endl; tree.in_orden(tree.root); cout << endl << endl;

	kd_node* target = new kd_node(point(d));
	vector<pair<kd_node*, double>> knn;

	tree.k_nearest_neighbor(target, k, tree.root, knn);
	cout << "target:             "; target->datum.printPoint();
	cout << "k nearest neighbor: " << endl;
	for (int i = 0; i < knn.size(); ++i)
	{
		knn[i].first->datum.printPoint(); cout << "distance: " << knn[i].second << endl;
	}

	cout << endl << "Distances from all nodes to the target: " << endl;

	for (int i = 0; i < size; ++i)
	{
		vec[i].printPoint();
		cout << vec[i].euclidean_distance(target->datum) << endl;
	}
}

void test_rq(int size, int d, bool random)
{
	kd_tree tree(d);
	vector<point> vec;

	cout << "Number of nodes: " << size << ", Dimension: " << d << endl << endl;
	cout << "Random nodes to insert:" << endl;
	for (int i = 0; i < size; ++i)
	{
		vec.push_back(point(d));
		vec.back().printPoint();
		tree.insert(vec.back());
	}
	cout << endl;

	cout << "In Orden:" << endl; tree.in_orden(tree.root); cout << endl;

	kd_node* target = new kd_node(point(d));
	vector<pair<kd_node*, double>> rq;

	int query_distance = 600;
	if (random) query_distance = generate_ramdon_number(500, 700);

	tree.range_query(target, query_distance, tree.root, rq);
	cout << "target:             "; target->datum.printPoint();
	cout << "query_distance: " << query_distance << endl;
	cout << "nodes results of query: " << endl;
	for (int i = 0; i < rq.size(); ++i)
	{
		rq[i].first->datum.printPoint(); cout << "distance: " << rq[i].second << endl;
	}

	cout << endl << "Distances from all nodes to the target: " << endl;

	for (int i = 0; i < size; ++i)
	{
		vec[i].printPoint();
		cout << vec[i].euclidean_distance(target->datum) << endl;
	}
}

int main()
{
	/*parameters(number od nodes, dimension of nodes)*/
	//test_insert(10, 3);

	/*parameters(number od nodes, dimension of nodes, number of nodes to erase)*/
	//test_erase(10, 4, 2);
	
	/*parameters(number od nodes, dimension of nodes, true = random point to serch or false = serch first random node that inserted)*/
	//test_serch(10, 3, true);
	
	/*parameters(number od nodes, dimension of nodes)*/
	//test_nn(10, 2);
	
	/*parameters(number od nodes, dimension of nodes, k value)*/
	test_knn(10, 2, 3);
	
	/*parameters(number od nodes, dimension of nodes, true = random query distance or false = predefine query distance)*/
	//test_rq(10, 3, false);
}