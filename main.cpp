#include "main.hpp"

using namespace Eigen;
using namespace std;


object * parse_test(istream &datafile);
vector<vector<float> > nums;
vector<vector<int> > facets;
vector<triangle> *triangles = new vector<triangle>();
vector<Vector3d> vertex_normal;

// given an index of vertex get the list of indices of
// triangles that have that vertex
unordered_map<int, vector<int> > tri_has_vertex;
//iterator for tri_has_vertex hash map
unordered_map<int, vector<int> >::iterator iter_vert;
int main(int argc, char** argv)
{	string tex_file;
	int  mesh_or_par;
	if (argc > 1) {
		mesh_or_par = atoi(argv[1]);
		cout<< mesh_or_par<<endl;
		tex_file = argv[2];
		cout<<tex_file;

	}	

    object *data = parse_test(cin);

    int num_face = (int) data->facets.size();

    nums = data->nums;
    facets = data->facets;

    GLfloat *vertex_coords = new GLfloat[num_face * 9];
    GLfloat *normal_coords = new GLfloat[num_face * 9];
    GLfloat *text_coords = new GLfloat[num_face * 9];

    for (int i = 0; i < num_face; i++) {

    	triangle t;
    	t.tri_num = i;
    	t.vertex_indices = vector<int> 
    		{facets[i][0], facets[i][1], facets[i][2]};


    	for (int j = 0; j < 3; j++) {
    		t.vertices_coords.push_back(nums[ facets[i][j] - 1 ]);
    		// try to find the vertex number and if it doesn't exist
    		if ((iter_vert = tri_has_vertex.find( facets[i][j] - 1) ) 
    			== tri_has_vertex.end() ) {
    			// make a vector
    			tri_has_vertex.insert( make_pair 
    				(facets[i][j] - 1, vector<int>{i}) );
    		}
    		// if it exists, just add it to the hash
    		else {
    			iter_vert->second.push_back(i);
    		}

    	}
    	t.calculate_normal();
    	//cout<<t.normal.transpose()<<endl;
    	triangles->push_back(t);
    }

    for (int i = 0; i < nums.size(); i++) {
    	Vector3d v;
    	v << 0, 0, 0;
    	vector<int> v2;
    	v2 = tri_has_vertex[i];
    	for (int j = 0; j < v2.size(); j++) {
    		v += triangles->at(v2[j]).normal;
    	}
    	v.normalize();

    	vertex_normal.push_back(v);
    }

    for (int i = 0; i < num_face; i++){
    	for (int j = 0; j < 3; j++) {
    		for (int k = 0; k < 3; k++) {
    		    vertex_coords[i * 9 + j * 3 + k] = 
    		    nums[facets[i][j] - 1][k];

    		    normal_coords[i * 9 + j * 3 + k] =
    		    vertex_normal[ facets[i][j] - 1](k);
    		}
    	}
    }

    /*cout<<tri_has_vertex.find(36)->second.size()<<endl;
    for (int i = 0; i < 3; i++) {
    	for(int j = 0; j < 3; j++){
	    	cout<<triangles->at(247).vertices_coords.at(i).at(j)<<endl;
	    }
    }*/
	
    std::fstream myfile(tex_file, std::ios_base::in);
 
    float a;
    int index = 0;
    while (myfile >> a)
    {
        text_coords[index++] = a;
        //if (index % 3 == 2) {
        //	text_coords[index++] = 0.0;
        //}
    }
    cout<<num_face*6<<endl;
    cout<<index<<endl;
    //assert(index == num_face *6);
    if (!mesh_or_par)
    	set_vertex(text_coords, num_face, normal_coords, text_coords,mesh_or_par);
    else
    	set_vertex (vertex_coords, num_face, normal_coords, text_coords,mesh_or_par);

    start(argc,argv);

    return 0;
}
