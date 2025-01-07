# 📈 Runge-Kutta 4

This project demonstrates the solution of **linear differential equations of any order** using the **4th-order Runge-Kutta method** (the 3/8 formula).

The main function, **`P1Z24_LSZ_runge-kutta.m`**, solves differential equations in the general form:
\[
a_m(x) \cdot y^{(m)} + a_{m-1}(x) \cdot y^{(m-1)} + \ldots + a_1(x) \cdot y' + a_0(x) \cdot y = b(x)
\]
by iteratively applying the 4th-order Runge-Kutta method.

---

## 🧪 Tests
- **`test*.m`**: Check the correctness of the implementation.
- **`test_num_*.m`**: Verify the numerical properties and performance of the method.

---

## 📄 Documentation
The project includes a detailed **report in Polish** (`P1Z24_LSZ.pdf`), which explains the theory behind the method, implementation details, and test results.
