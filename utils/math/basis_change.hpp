#include "utils/Vector.hpp"

#include <tuple>
namespace Utils {
namespace Math {
inline 
Vector3d basis_change(const std::tuple<Vector3d, Vector3d, Vector3d>& basis_vectors,
                      const Vector3d& v) {
    // Transformation matrix has new basis vectors as columns
    Vector3d row_0 = Vector3d{{std::get<0>(basis_vectors)[0],std::get<1>(basis_vectors)[0], std::get<2>(basis_vectors)[0]}};
    Vector3d row_1 = Vector3d{{std::get<0>(basis_vectors)[1],std::get<1>(basis_vectors)[1], std::get<2>(basis_vectors)[1]}};
    Vector3d row_2 = Vector3d{{std::get<0>(basis_vectors)[2],std::get<1>(basis_vectors)[2], std::get<2>(basis_vectors)[2]}};

    return Vector3d{{row_0* v, row_1 * v, row_2 * v}};
}

}
}
