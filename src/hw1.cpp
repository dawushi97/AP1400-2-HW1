#include "hw1.h"
#include <random>
#include <vector>

// user-defined functions
int algebra::rowNum(const Matrix &matrix) {
  if (matrix.empty()) {
    return 0;
  }
  return matrix.size();
}

int algebra::colNum(const Matrix &matrix) {
  if (matrix.empty()) {
    return 0;
  }
  return matrix[0].size();
}

bool algebra::isEmpty(const Matrix &matrix) { return matrix.empty(); }

double &algebra::element(Matrix &matrix, int row, int col) {
  return matrix[row][col];
}

double algebra::dotProduct(const Matrix &matrix1, const Matrix &matrix2,
                           int row, int col) {
  double dotProduct = 0;
  for (size_t i = 0; i < matrix1.size(); i++) {
    dotProduct += matrix1[row][i] * matrix2[i][col];
  }
  return dotProduct;
}

Matrix algebra::adjugate(const Matrix &matrix) {
  Matrix result(matrix.size(), std::vector<double>(matrix[0].size(), 0));

  // Calculate cofactor matrix
  for (size_t i = 0; i < matrix.size(); i++) {
    for (size_t j = 0; j < matrix[0].size(); j++) {
      // Add sign term (-1)^(i+j)
      double sign = ((i + j) % 2 == 0) ? 1 : -1;
      result[i][j] = sign * determinant(minor(matrix, i, j));
    }
  }
  return algebra::transpose(result);
}
//Hw Need to implement
Matrix algebra::zeros(size_t n, size_t m) {
  Matrix matrix(n, std::vector<double>(m, 0));
  return matrix;
}
Matrix algebra::ones(size_t n, size_t m) {
  Matrix matrix(n, std::vector<double>(m, 1));
  return matrix;
}
Matrix algebra::random(size_t n, size_t m, double min, double max) {
  if (min > max) {
    throw logic_error("min > max");
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(min, max);

  Matrix matrix(n, std::vector<double>(m));
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; j++) {
      matrix[i][j] = dis(gen);
    }
  }
  return matrix;
}
void algebra::show(const Matrix &matrix) {

  cout << std::fixed;           // 使用固定小数点
  cout << std::setprecision(3); // 设置小数点后位数

  for (const auto &row : matrix) {
    for (double val : row) {
      // setw(8): 设置字段宽度为8
      // right: 右对齐
      cout << std::setw(8) << std::right << val << " ";
    }
    cout << endl;
  }
}

Matrix algebra::multiply(const Matrix &matrix, double c) {
  size_t n = matrix.size();
  size_t m = matrix[0].size();
  Matrix result(n, std::vector<double>(m, 0));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; ++j) {
      result[i][j] = matrix[i][j] * c;
    }
  }
  return result;
}

Matrix algebra::multiply(const Matrix &matrix1, const Matrix &matrix2) {
  size_t cow1 = matrix1.size();
  size_t cow2 = matrix2.size();
  size_t com1 = matrix1[0].size();
  size_t com2 = matrix2[0].size();
  if (com1 != cow2) {
    Matrix result(cow1, std::vector<double>(com2, 0));
    for (size_t i = 0; i < cow1; i++) {
      for (size_t j = 0; j < cow2; j++) {
        result[i][j] = dotProduct(matrix1, matrix2, i, j);
      }
    }
  }
}

Matrix algebra::sum(const Matrix &matrix, double c) {
  Matrix result(matrix.size(), std::vector<double>(matrix[0].size(), 0));
  for (size_t i = 0; i < matrix.size(); i++) {
    for (size_t j = 0; j < matrix[0].size(); j++) {
      result[i][j] = matrix[i][j] + c;
    }
  }
  return result;
}

Matrix algebra::sum(const Matrix &matrix1, const Matrix &matrix2) {
  Matrix result(matrix1.size(), std::vector<double>(matrix1.size(), 0));
  for (size_t i = 0; i < matrix1.size(); i++) {
    for (size_t j = 0; j < matrix1.size(); j++) {
      result[i][j] = matrix1[i][j] + matrix2[i][j];
    }
  }
  return result;
}

Matrix algebra::transpose(const Matrix &matrix) {
  int result_col = rowNum(matrix);
  int result_row = colNum(matrix);
  Matrix result(result_row, std::vector<double>(result_col, 0));
  for (int i = 0; i < result_col; i++) {
    for (int j = 0; j < result_row; j++) {
      result[j][i] = matrix[i][j];
    }
  }
  return result;
}

Matrix algebra::minor(const Matrix &matrix, size_t n, size_t m) {
  int result_row = rowNum(matrix) - 1;
  int result_col = colNum(matrix) - 1;
  Matrix result(result_row, std::vector<double>(result_col, 0));

  int new_i = 0;
  for (size_t i = 0; i < rowNum(matrix); i++) {
    if (i == n)
      continue;
    int new_j = 0;
    for (size_t j = 0; j < colNum(matrix); j++) {
      if (j == m)
        continue;
      result[new_i][new_j] = matrix[i][j];
      new_j++;
    }
    new_i++;
  }
  return result;
}

double algebra::determinant(const Matrix &matrix) {
  int n = rowNum(matrix);
  if (n == 1) {
    return matrix[0][0];
  } else if (n == 2) {
    return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
  }
  double det = 0;
  for (int i = 0; i < n; i++) {
    double cofactor = ((i % 2 == 0) ? 1 : -1) * matrix[0][i] *
                      determinant(minor(matrix, i, n));
    det += cofactor;
  }
  return det;
}

Matrix algebra::inverse(const Matrix &matrix) {
  if (isEmpty(matrix)) {
    return matrix;
  }
  if (colNum(matrix) != rowNum(matrix)) {
    throw logic_error("Matrix is not square");
  }
  if (determinant(matrix) == 0) {
    throw logic_error("Matrix is not invertible");
  }
  return multiply(adjugate(matrix), 1 / determinant(matrix));
}

Matrix algebra::concatenate(const Matrix &matrix1, const Matrix &matrix2,
                            int axis) {
  if (axis == 0) {
    if (colNum(matrix1) != rowNum(matrix2)) {
      throw logic_error("Matrix is not square");
    }
    int total_rows = rowNum(matrix1) + rowNum(matrix2);
    int col = colNum(matrix1);
    Matrix result(total_rows, std::vector<double>(col, 0));
    for (int i = 0; i < rowNum(matrix1); i++) {
      for (int j = 0; j < col; j++) {
        result[i][j] = matrix1[i][j];
      }
    }
    for (int i = 0; i < rowNum(matrix2); i++) {
      for (int j = 0; j < col; j++) {
        result[i + rowNum(matrix1)][j] = matrix2[i][j];
      }
    }
    return result;
  } else if (axis == 1) {
    if (rowNum(matrix1) != colNum(matrix2)) {
      throw logic_error("Matrix is not square");
    }
    int rows = rowNum(matrix1);
    int total_cols = colNum(matrix2) + colNum(matrix1);
    Matrix result(rows, std::vector<double>(total_cols, 0));
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < colNum(matrix1); j++) {
        result[i][j] = matrix1[i][j];
      }
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < colNum(matrix2); j++) {
        result[i][j + colNum(matrix1)] = matrix2[i][j];
      }
    }
    return result;
  }
}

Matrix algebra::ero_swap(const Matrix &matrix, size_t r1, size_t r2) {
  if (r1 > rowNum(matrix) || r2 > colNum(matrix)) {
    throw logic_error("Matrix is not square");
  }
  Matrix result = zeros(rowNum(matrix), colNum(matrix));
  std::swap(result[r1], result[r2]);
  return result;
}

Matrix algebra::ero_multiply(const Matrix &matrix, size_t r, double c) {
  if (r >= rowNum(matrix)) {
    throw logic_error("Matrix is not square");
  }
  Matrix result = zeros(rowNum(matrix), colNum(matrix));
  for (size_t j = 0; j < rowNum(matrix); j++) {
    result[r][j] *= c;
  }
  return result;
}

Matrix algebra::ero_sum(const Matrix &matrix, size_t r1, double c, size_t r2) {
  if (r1 > rowNum(matrix) || r2 > rowNum(matrix)) {
    throw logic_error("Matrix is not square");
  }
  Matrix result = zeros(rowNum(matrix), colNum(matrix));
  for (size_t j = 0; j < rowNum(matrix); j++) {
    result[r2][j] += result[r1][j] * c;
  }
  return result;
}

Matrix algebra::upper_triangular(const Matrix &matrix) {
  Matrix result = matrix; // Create a copy to modify
  size_t n = matrix.size();

  // 对每一列进行操作
  for (size_t i = 0; i < n - 1; i++) {
    // 找到当前列中最大元素的行（部分主元消去法）
    size_t max_row = i;
    for (size_t k = i + 1; k < n; k++) {
      if (std::abs(result[k][i]) > std::abs(result[max_row][i])) {
        max_row = k;
      }
    }

    // 如果需要，交换行
    if (max_row != i) {
      std::swap(result[i], result[max_row]);
    }

    // 如果主元为0，继续下一列
    if (std::abs(result[i][i]) < 1e-10) {
      continue;
    }

    // 消元操作
    for (size_t j = i + 1; j < n; j++) {
      double factor = result[j][i] / result[i][i];
      for (size_t k = i; k < n; k++) {
        result[j][k] -= factor * result[i][k];
      }
    }
  }
  return result;
}

