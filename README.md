# WinQC - VQE on DAQC

*A challenge provided by IQM Quantum Computers*

This repository contains the WinQC team's solution to the DAQC VQE challenge for the Womanium Hackathon 2022.

**Team Name**: <font color='green'>**WinQC**</font>

**Team members**:

Nikolas Klemola Tango (Discord id-> nick#2491 , Github id-> Nickname88888, Email-> nikoklemola2011@hotmail.com)

Tereza Viskova (Discord id-> Terkasaurus#6388 , Github id-> tervisk, Email-> terezaviskova1997@gmail.com )

Mohammad Zuhair Khan (Discord id-> Raptor#6814 , Github id-> MZuhairKhan, Email-> khanmohammadzuhair@gmail.com)

César Bertoni Ocampo (Discord id-> Ceboc#2938 , Github id->Ceboc, Email-> bertoni@ciencias.unam.mx) 

Shilan Abo   (Discord id-> shilqc #7311, Github id-> Shilqc, Email-> shsavan@icloud.com) 

**Pitch Presenter**: Mohammad Zuhair Khan

**Challenge:**: Digital-analog Variational Quantum Eigensolver (IQM)

**Overview**

The Variational Quantum Eigensolver algorithm (VQE) method provides a way of approximating the lowest energy eigenstate. In this approach we prepare a guess state $\psi(\Theta)$ as $|\psi(\Theta)\rangle = \sum_n a_n |\psi_n\rangle$.

Then, for any arbitrary state, we know 

$$\langle \psi(\Theta)|H |\psi(\theta)\rangle = (\sum_n a^*_n \langle \psi_n|)H(\sum_m a_m |\psi_m\rangle) = \sum_n |a_n|^2 \langle \psi_n|H |\psi_n\rangle= \sum a_n^2 E_n \ge \sum_n a_n^2 E_0 \approx E_0$$

So in Variational principle $\frac{\langle \psi(\Theta)| H |\psi(\Theta)\rangle}{\langle \psi(\Theta)|\psi(\Theta)\rangle} \ge E_0$, where here $E_{guess}=\langle \psi(\Theta)| H |\psi(\Theta)\rangle$.

We can optimize $\Theta$ in the guess state so that we are left (hopefully) with a close upper bound to $E_0$.

**Implementation**

Because gate noise is a problem in the NISQ era, we combine digital single-qubit operations with analog multi-qubit entangling blocks in a method known as digital-analog quantum computing (DAQC). We can use either a homogeneous nearest neighbor Hamiltonian $(H_{zz}=\sum_{j}g_{jj+1}Z^jZ^{j+1})$ or a homogeneous all-to-all two-body Ising Hamiltonian $(H_{zz}=\sum_{j>k}g_{jk}Z^jZ^k)$ as analog blocks.

1. The first step is to map the fermionic Hamiltonian of a given molecule into the qubit Hamiltonian.
2. Then we need to prepare the guess state, to do that we create/update an ansatz for state preparation on the quantum computer and build the quantum circuit by combining the entangling operations under a given Hamiltonian with single qubit gates. 
3. Then, on the basis of our qubit Hamiltonian, calculate the state's expectation values.
4. Finally send the results to a classical optimizer to update the gate/wave parameters and repeat steps 2-4 until convergence.

In our approach we study Hydrogen molecule $H_2$ in a two-qubit case.

In a VQE method, the quantum circuit used for a subroutine is called ansatz. These can consist of single qubit digital gates or entanglement analog gates. In the DAQC VQE method, an analog block ansatz is used in the form of Ising model Hamiltonian.

We used Simultaneous Perturbation Stochastic Approximation (SPSA) optimizer, which is ideal for “noisy” cost functions. 

**Effect of DAQC**

The choice of Hamiltonian affects the degrees of freedom and the SQRs. 

The benefit of DAQC VQE compared to a method with only digital blocks is the ability of entanglement in the analog ansatz. This provides significant speed up and more iterations, thus more accurate result. 

### **Results**
The tests included three different types of DAQC hamiltonian implementations along with a purely digital implementation of a VQE as a comparison.

The three different 2-qubit Hamiltonian implementation types were the following:
 1. Both time and the hamiltonian coefficients were parameters, but the coefficients were constant between the analog blocks.
 2. Both time and the hamiltonian coefficients were parameters, and the coefficients could vary between the analog blocks.
 3. The hamiltonian coefficients were fixed as $1.0$ and time was the only parameter.

A purely digital implementation of VQE is also used as a comparison.

Below are the graphs comparing the three different implementations and how well they are able to estimate the ground state value, the x-axis shows how these results change when more digital-analog blocks are added. The digital-only method of estimation is absent from the graph as it is several orders of magnitude further away from the reference value at $\Delta = 0.79027$.

![image.png](https://i.imgur.com/jFRtDwe.png)

