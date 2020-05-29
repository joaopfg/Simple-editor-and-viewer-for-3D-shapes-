#include <igl/opengl/glfw/Viewer.h>
#include <ostream>
#include <math.h>

using namespace Eigen;
using namespace std;

/**
 * A class for representing linear transformations on 3D points (using homogeneous coordinates)
 * */
class MeshTransformation
{
public:
/*
Initialize the identity transformation
**/
  MeshTransformation()
  {
    MatrixXd m(4, 4);
    m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
    m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    M = m;
  }


/*
Initialize a scaling transformation
**/
  MeshTransformation(double s1, double s2, double s3)
  {
    //COMPLETED

    MatrixXd m(4, 4);
    m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
    m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    m.row(0) *= s1;
    m.row(1) *= s2;
    m.row(2) *= s3;

    M = m;
  }

/*
Initialize a rotation transformation around a given axis (X, Y or Z) <br><br>

 @param  direction  a value 0, 1 or 2 indicating the direction (X, Y or Z respectively)
**/
  MeshTransformation(double theta, int direction)
  {
    //COMPLETED

    MatrixXd m(4, 4);
    m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
    m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    if(direction == 0){
      m(1, 1) = cos(theta);
      m(1, 2) = -sin(theta);
      m(2, 1) = sin(theta);
      m(2, 2) = cos(theta);
    }
    else if(direction == 1){
      m(0, 0) = cos(theta);
      m(0, 2) = sin(theta);
      m(2, 0) = -sin(theta);
      m(2, 2) = cos(theta);
    }
    else{
      m(0, 0) = cos(theta);
      m(0, 1) = -sin(theta);
      m(1, 0) = sin(theta);
      m(1, 1) = cos(theta);
    }

    M = m;
  }

/*
Initialize a translation
**/
  MeshTransformation(RowVector3d t)
  {
    //COMPLETED

    MatrixXd m(4, 4);
    m(0, 0) = 1.0; m(1, 0) = 0.0; m(2, 0) = 0.0; m(3, 0) = 0.0;
    m(0, 1) = 0.0; m(1, 1) = 1.0; m(2, 1) = 0.0; m(3, 1) = 0.0;
    m(0, 2) = 0.0; m(1, 2) = 0.0; m(2, 2) = 1.0; m(3, 2) = 0.0;
    m(0, 3) = 0.0; m(1, 3) = 0.0; m(2, 3) = 0.0; m(3, 3) = 1.0;

    m(0, 3) = t(0, 0);
    m(1, 3) = t(0, 1);
    m(2, 3) = t(0, 2);

    M = m;
  }

/*
Matrix accessor

@return  the matrix transformation
**/
  MatrixXd get_matrix() {
    return M;
  }

/*
Initialize a transformation given an input matrix 'm'
**/
  void set_matrix(MatrixXd m)
  {
    M = m;
  }

/*
Apply the transformation to all vertices stored in a matrix 'V' <br>

@param V  vector storing the input points
**/
  void transform(MatrixXd &V) {
    //COMPLETED
    int n = V.rows();

    for(int i=0;i<n;i++) V.row(i) = transform(V.row(i));
  }

  	/**
	 * Apply the transformation to a 3d (row) vector 'v' <br>
   * 
   * Remark: use homogeneous coordinates
   * 
   * @return  the vector after transformation
	 */
	RowVector3d transform(RowVector3d v) {
    //COMPLETED
    RowVector3d result;
    Vector4d Homogeneous_coord(v(0, 0), v(0, 1), v(0, 2), 1.0);

    Homogeneous_coord = M*Homogeneous_coord;

    result(0, 0) = Homogeneous_coord(0, 0);
    result(0, 1) = Homogeneous_coord(1, 0);
    result(0, 2) = Homogeneous_coord(2, 0);

		return result;
	}

	/**
	 * Compose the current transformation with a transfomation 't': return a new transformation
	 */
	MeshTransformation compose(MeshTransformation t) {
    //COMPLETED
    MeshTransformation res(1.0, 1.0, 1.0);
    MatrixXd m = res.get_matrix();
    m = (t.get_matrix())*m;
    res.set_matrix(m);

    return res;
	}

	/**
	 * Print the matrix transformation
	 */
  friend std::ostream& operator<<(std::ostream &os, MeshTransformation& t) {
    return os << "matrix:\n" << t.get_matrix() << std::endl;
  }

private:
  MatrixXd M; // a 4x4 matrix representing a linear transformation
};
