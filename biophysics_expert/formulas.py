"""
Formula definitions for the BiophysicsExpert system.

This module contains the default formulas included with the package.
Users can extend this list with custom formulas.
"""

# Default formulas included in the package
DEFAULT_FORMULAS = [
    {
                "Formula_Name": "First Law of Thermodynamics",
                "Unknown": "dE",
                "Formula": "dE = δQ + δW",
                "LaTeX": "dE = \\delta Q + \\delta W",
                "Variables": ["δQ", "δW"],
                "Additional_Info": "Where δ denotes an infinitesimal amount, Q is heat, and W is work."
            },
            {
                "Formula_Name": "First Law of Thermodynamics for Pressure-Volume Work",
                "Unknown": "dE",
                "Formula": "dE = δQ - pdV",
                "LaTeX": "dE = \\delta Q - p\\,dV",
                "Variables": ["δQ", "p", "dV"],
                "Additional_Info": "Where δ denotes an infinitesimal amount, Q is heat, p is pressure, and dV is change in volume."
            },
            {
                "Formula_Name": "First Law of Thermodynamics for Reversible Processes",
                "Unknown": "dE",
                "Formula": "dE = TdS - pdV",
                "LaTeX": "dE = T\\,dS - p\\,dV",
                "Variables": ["T", "dS", "p", "dV"],
                "Additional_Info": "Where T is temperature, p is pressure, dS is change in entropy, and dV is change in volume."
            },
            {
                "Formula_Name": "First Law of Thermodynamics (Spontaneous Processes and Isolated Systems)",
                "Unknown": "ΔE",
                "Formula": "ΔE = E_{eq} - E_{neq}",
                "LaTeX": "\\Delta E = E_{eq} - E_{neq}",
                "Variables": ["E_{eq}", "E_{neq}"],
                "Additional_Info": "where ΔE < 0 for spontaneous processes and dE = 0 for isolated systems."
            },

            # Entropy Differential Equation
            {
                "Formula_Name": "Entropy Differential Equation",
                "Unknown": "dS",
                "Formula": "dS = (∂S/∂E)_V dE + (∂S/∂V)_E dV",
                "LaTeX": "dS = \\left(\\frac{\\partial S}{\\partial E}\\right)_V dE + \\left(\\frac{\\partial S}{\\partial V}\\right)_E dV",
                "Variables": ["∂S", "∂E", "V", "dE", "∂V", "E", "dV"],
                "Additional_Info": "Where ∂ is partial derivative, S is entropy, E is energy, V is volume, dE is change in energy, and dV is change in volume."
            },
            {
                "Formula_Name": "Entropy Differential Equation",
                "Unknown": "dS",
                "Formula": "dS = (1/T)dE + (∂S/∂V)_E dV",
                "LaTeX": "dS = \\frac{1}{T}dE + \\left(\\frac{\\partial S}{\\partial V}\\right)_E dV",
                "Variables": ["T", "dE", "∂S", "∂V", "E", "dV"],
                "Additional_Info": "Where ∂ is a partial derivative, T is temperature, S is entropy, V is volume, and dE and dV are changes in energy and volume."
            },
            {
                "Formula_Name": "Entropy Differential Equation For Reversible Processes",
                "Unknown": "dS",
                "Formula": "dS = δQ/T",
                "LaTeX": "dS = \\frac{\\delta Q}{T}",
                "Variables": ["δQ", "T"],
                "Additional_Info": "Where δ denotes an infinitesimal amount, Q is heat, and T is temperature."
            },
            {
                "Formula_Name": "Entropy Differential Equation (Expanded Form)",
                "Unknown": "dS",
                "Formula": "dS = δQ/T + [-P/T + (∂S/∂V)_E]dV",
                "LaTeX": "dS = \\frac{\\delta Q}{T} + \\left[-\\frac{P}{T} + \\left(\\frac{\\partial S}{\\partial V}\\right)_E\\right]dV",
                "Variables": ["δQ", "T", "P", "∂S", "∂V", "E", "dV"],
                "Additional_Info": "Where δ denotes an infinitesimal amount, ∂ is a partial derivative, Q is heat, T is temperature, p is pressure, S is entropy, V is volume, E is energy, and dV is change in volume."
            },
            {
                "Formula_Name": "Entropy Differential Equation (Spontaneous Processes and at Equilibrium)",
                "Unknown": "ΔS",
                "Formula": "ΔS = S_{eq} - S_{neq}",
                "LaTeX": "\\Delta S = S_{eq} - S_{neq}",
                "Variables": ["S_{eq}", "S_{neq}"],
                "Additional_Info": "Where ΔS>0 for spontaneous processes and ΔS=0 at equilibrium."
            },

            # Fundamental Gibbs Equation
            {
                "Formula_Name": "Fundamental Gibbs Equation",
                "Unknown": "dE",
                "Formula": "dE = TdS - pdV + μdN",
                "LaTeX": "dE = T\\,dS - p\\,dV + \\mu\\,dN",
                "Variables": ["T", "dS", "p", "dV", "μ", "dN"],
                "Additional_Info": "Where T is temperaute, p is pressure, μ is chemical potential, and dS, dV, and dN are the changes in entropy, volume, and the number of particles."
            },
            {
                "Formula_Name": "Fundamental Gibbs Equation In Terms of Enthalpy",
                "Unknown": "dH",
                "Formula": "dH = TdS + Vdp + μdN",
                "LaTeX": "dH = T\\,dS + V\\,dp + \\mu\\,dN",
                "Variables": ["T", "dS", "V", "dp", "μ", "dN"],
                "Additional_Info": "Where T is temperature, V is volume, μ is chemical potential, and dS, dP, and dN are the changes in entropy, pressure, and the number of particles."
            },
            {
                "Formula_Name": "Fundamental Gibbs Equation In Terms of Gibbs Free Energy",
                "Unknown": "dG",
                "Formula": "dG = -SdT + Vdp + μdN",
                "LaTeX": "dG = -S\\,dT + V\\,dp + \\mu\\,dN",
                "Variables": ["S", "dT", "V", "dp", "μ", "dN"],
                "Additional_Info": "Where S is entropy, V is volume, μ is chemical potential, and dT, dp, and dN are the changes in temperature, pressure, and the number of particles."
            },
            {
                "Formula_Name": "Fundamental Gibbs Equation In Terms of Helmholtz Free Energy",
                "Unknown": "dA",
                "Formula": "dA = -SdT - pdV + μdN",
                "LaTeX": "dA = -S\\,dT - p\\,dV + \\mu\\,dN",
                "Variables": ["S", "dT", "p", "dV", "μ", "dN"],
                "Additional_Info": "Where S is entropy, p is pressure, μ is chemical potential, and dT, dV, and dN are the changes in temperature, volume, and the number of particles."
            },

            # Heat Transfer and Entropy Change
            {
                "Formula_Name": "Heat Transfer and Entropy Change",
                "Unknown": "dS",
                "Formula": "dS = δQ/T",
                "LaTeX": "dS = \\frac{\\delta Q}{T}",
                "Variables": ["δQ", "T"],
                "Additional_Info": "Where δ denotes an infinitesimal amount, Q is heat, and T is temperature."
            },
            {
                "Formula_Name": "Heat Transfer and Entropy Change For Variable Temperature Processes",
                "Unknown": "ΔS",
                "Formula": "ΔS = ∫δQ/T",
                "LaTeX": "\\Delta S = \\int\\frac{\\delta Q}{T}",
                "Variables": ["∫", "δQ", "T"],
                "Additional_Info": "Where ∫ denotes an intergral, δ denotes an infinitesimal amount, Q is heat, and T is temperature."
            },
            {
                "Formula_Name": "Heat Transfer and Entropy Change For Constant Temperature Processes",
                "Unknown": "ΔS",
                "Formula": "ΔS = Q/T",
                "LaTeX": "\\Delta S = \\frac{Q}{T}",
                "Variables": ["Q", "T"],
                "Additional_Info": "Where Q is heat and T is temperature."
            },

            # Entropy Variation in Connected Subsystems
            {
                "Formula_Name": "Entropy Variation in Connected Subsystems",
                "Unknown": "δS",
                "Formula": "δS = δS_1 + δS_2",
                "LaTeX": "\\delta S = \\delta S_1 + \\delta S_2",
                "Variables": ["δS1", "δS2"],
                "Additional_Info": "Where δ denotes an infinitesimal amount, and S1 and S2 are the entropies in connected subsystems."
            },
            {
                "Formula_Name": "Entropy Variation in Connected Subsystems",
                "Unknown": "δS",
                "Formula": "δS = (∂S_1/∂E_1)δE_1 + (∂S_2/∂E_2)δE_2",
                "LaTeX": "\\delta S = \\left(\\frac{\\partial S_1}{\\partial E_1}\\right)\\delta E_1 + \\left(\\frac{\\partial S_2}{\\partial E_2}\\right)\\delta E_2",
                "Variables": ["∂S1", "∂E1", "δE1", "∂S2", "∂E2", "δE2"],
                "Additional_Info": "Where δ denotes an infinitesimal amount, ∂ is a partial derivative, S1 and S2 are the entropies in connected subsystems, and E1 and E2 are energies in connected subsystems."
            },
            {
                "Formula_Name": "Entropy Variation in Connected Subsystems",
                "Unknown": "δS",
                "Formula": "δS = (1/T_1)δE_1 + (1/T_2)δE_2",
                "LaTeX": "\\delta S = \\frac{1}{T_1}\\delta E_1 + \\frac{1}{T_2}\\delta E_2",
                "Variables": ["T1", "δE1", "T2", "δE2"],
                "Additional_Info": "Where T1 and T2 are temperatures in connected subsystems and E1 and E2 are energies in connected subsystems."
            },
            {
                "Formula_Name": "Entropy Variation in Connected Subsystems",
                "Unknown": "δS",
                "Formula": "δS = (1/T_1 - 1/T_2)δE_1",
                "LaTeX": "\\delta S = \\left(\\frac{1}{T_1}-\\frac{1}{T_2}\\right)\\delta E_1",
                "Variables": ["T1", "T2", "δE1"],
                "Additional_Info": "Where T1 and T2 are temperatures in connected subsystems and δE1 = -δE2."
            },

            # Helmholtz Free Energy
            {
                "Formula_Name": "Helmholtz Free Energy",
                "Unknown": "F",
                "Formula": "F = E - TS",
                "LaTeX": "F = E - TS",
                "Variables": ["E", "T", "S"],
                "Additional_Info": "Where E is energy, T is temperature, and S is entropy."
            },
            {
                "Formula_Name": "Helmholtz Free Energy",
                "Unknown": "F",
                "Formula": "F = U - TS",
                "LaTeX": "F = U - TS",
                "Variables": ["U", "T", "S"],
                "Additional_Info": "Where U is internal energy, T is temperature, and S is entropy."
            },
            {
                "Formula_Name": "Helmholtz Free Energy (Incorporating State Variables)",
                "Unknown": "F",
                "Formula": "F = -pV + μN - TS",
                "LaTeX": "F = -pV + \\mu N - TS",
                "Variables": ["p", "V", "μ", "N", "T", "S"],
                "Additional_Info": "Where p is pressure, V is volume, μ is chemical potential, N is the number of particles in the system, T is temperature, and S is entropy."
            },
            {
                "Formula_Name": "RNA Free Energy Function",
                "Unknown": "F",
                "Formula": "F = -(1/β)lnZ",
                "LaTeX": "F = -\\frac{1}{\\beta}\\ln Z",
                "Variables": ["β", "Z"],
                "Additional_Info": "Where β = 1/kBT (kB is the Boltzman constant) and Z is the partition function."
            },
            {
                "Formula_Name": "Helmholtz Free Energy",
                "Unknown": "F",
                "Formula": "F = -k_BTlnZ",
                "LaTeX": "F = -k_B T \\ln Z",
                "Variables": ["kB", "T", "Z"],
                "Additional_Info": "Where kB is the Boltzman constant, T is temperature, and Z is the partition function."
            },

            # Helmholtz Free Energy Change
            {
                "Formula_Name": "Helmholtz Free Energy Change",
                "Unknown": "dF",
                "Formula": "dF = dE - TdS - SdT",
                "LaTeX": "dF = dE - T\\,dS - S\\,dT",
                "Variables": ["dE", "T", "dS", "S", "dT"],
                "Additional_Info": "Where T is temperature, S is entropy, and dE, dS, and dT are the changes in energy, entropy, and temperature respectively."
            },
            {
                "Formula_Name": "Helmholtz Free Energy Change",
                "Unknown": "ΔF",
                "Formula": "ΔF = ΔE - TΔS",
                "LaTeX": "\\Delta F = \\Delta E - T\\Delta S",
                "Variables": ["ΔE", "T", "ΔS"],
                "Additional_Info": "Where T is temperature and ΔE and ΔS are the changes in energy and entropy respectively."
            },
            {
                "Formula_Name": "Helmholtz Free Energy Change (Natural Variable Form)",
                "Unknown": "dF",
                "Formula": "dF = -SdT - pdV + μdN",
                "LaTeX": "dF = -S\\,dT - p\\,dV + \\mu\\,dN",
                "Variables": ["S", "dT", "p", "dV", "μ", "dN"],
                "Additional_Info": "Where S is entropy, p is pressure, μ is chemical potential, and dT, dV, and dN are the changes in temperature, volume, and the number of particles."
            },
            {
                "Formula_Name": "Helmholtz Free Energy Change (For Closed Systems)",
                "Unknown": "dF",
                "Formula": "dF = -SdT - pdV",
                "LaTeX": "dF = -S\\,dT - p\\,dV",
                "Variables": ["S", "dT", "p", "dV"],
                "Additional_Info": "Where S is entropy, p is pressure, dT is the change in temperature and dV is the change in volume."
            },

            # Gibbs Free Energy
            {
                "Formula_Name": "Gibbs Free Energy",
                "Unknown": "G",
                "Formula": "G = F + pV",
                "LaTeX": "G = F + pV",
                "Variables": ["F", "p", "V"],
                "Additional_Info": "Where F is the Helmholtz free energy and pV."
            },
            {
                "Formula_Name": "Gibbs Free Energy",
                "Unknown": "G",
                "Formula": "G = E - TS + pV",
                "LaTeX": "G = E - TS + pV",
                "Variables": ["E", "T", "S", "p", "V"],
                "Additional_Info": "Where E is energy, T is temperature, S is entropy, and pV is the change in volume."
            },
            {
                "Formula_Name": "Gibbs Free Energy",
                "Unknown": "G",
                "Formula": "G = H - TS",
                "LaTeX": "G = H - TS",
                "Variables": ["H", "T", "S"],
                "Additional_Info": "Where H is enthalpy (defined as E+pV), T is temperature, and S is entropy."
            },

            # Gibbs Free Energy Change
            {
                "Formula_Name": "Gibbs Free Energy Change",
                "Unknown": "dG",
                "Formula": "dG = dF + pdV + Vdp",
                "LaTeX": "dG = dF + p\\,dV + V\\,dp",
                "Variables": ["dF", "p", "dV", "V", "dp"],
                "Additional_Info": "Where p is pressure, dV is change in volume, and dF is change in Helmholtz free energy."
            },
            {
                "Formula_Name": "Gibbs Free Energy Change (Natural Variables Form)",
                "Unknown": "dG",
                "Formula": "dG = -SdT + Vdp + μdN",
                "LaTeX": "dG = -S\\,dT + V\\,dp + \\mu\\,dN",
                "Variables": ["S", "dT", "V", "dp", "μ", "dN"],
                "Additional_Info": "Where S is entropy, V is volume, μ is chemical potential, and dT, dp, and dN are the changes in temperature, pressure, and the number of particles."
            },
            {
                "Formula_Name": "Gibbs Free Energy Change (For Closed Systems)",
                "Unknown": "dG",
                "Formula": "dG = -SdT + Vdp",
                "LaTeX": "dG = -S\\,dT + V\\,dp",
                "Variables": ["S", "dT", "V", "dp"],
                "Additional_Info": "Where S is entropy, V is volume, dT is the change in temperature, and dp is the change in pressure."
            },

            # Microcanonical Ensemble Probability
            {
                "Formula_Name": "Microcanonical Ensemble Probability",
                "Unknown": "Pn",
                "Formula": "P_n = 1/Ω(N,V,E)",
                "LaTeX": "P_n = \\frac{1}{\\Omega(N,V,E)}",
                "Variables": ["Ω", "N", "V", "E"],
                "Additional_Info": "Where Ω is the total number of microstates, N is the number of molecules, V is volume, and E is the energy."
            },

            # Boltzmann Probability Distribution
            {
                "Formula_Name": "Boltzmann Probability Distribution",
                "Unknown": "Pn",
                "Formula": "P_n = e^(-βE_n)/Z",
                "LaTeX": "P_n = \\frac{e^{-\\beta E_n}}{Z}",
                "Variables": ["β", "En", "Z"],
                "Additional_Info": "Where β is equal to 1/kBT, En is the energy of state n, and Z is the partition function."
            },
            {
                "Formula_Name": "Boltzmann Probability Distribution (Using Helmholtz Free Energy)",
                "Unknown": "Pn",
                "Formula": "P_n = e^(β(F-E_n))",
                "LaTeX": "P_n = e^{\\beta(F-E_n)}",
                "Variables": ["β", "F", "En"],
                "Additional_Info": "Where β is equal to 1/kBT, F is the Helmholtz free energy, and En is the energy of state n."
            },
            {
                "Formula_Name": "Boltzmann Probability Distribution (Derived from Micocanonical Ensemble)",
                "Unknown": "Pn",
                "Formula": "P_n = Ω(E_B+E-E_n)/Ω(E_0)",
                "LaTeX": "P_n = \\frac{\\Omega(E_B+E-E_n)}{\\Omega(E_0)}",
                "Variables": ["Ω", "EB", "E", "En", "E0"],
                "Additional_Info": "Where Ω is the number of microstates or density of states, EB is the energy of a bath of reservoir, En is the energy of state n, and E0 is the reference energy (such as the total energy of the system)."
            },

            # Statistical Definition of Entropy
            {
                "Formula_Name": "Statistical Definition of Entropy",
                "Unknown": "S",
                "Formula": "S = k_BlnΩ(N,V,E)",
                "LaTeX": "S = k_B\\ln\\Omega(N,V,E)",
                "Variables": ["kB", "Ω", "N", "V", "E"],
                "Additional_Info": "Where kB is the Boltzman constant, Ω is the total number of microstates, N is the number of molecules, V is volume, and E is the energy."
            },
            {
                "Formula_Name": "Entropy",
                "Unknown": "S",
                "Formula": "S = k_Blng(E)",
                "LaTeX": "S = k_B\\ln g(E)",
                "Variables": ["kB", "g(E)"],
                "Additional_Info": "Where kB is the Boltzman constant and g(E) is the density of states."
            },
            {
                "Formula_Name": "Entropy",
                "Unknown": "S",
                "Formula": "S/k_B = lnΩ",
                "LaTeX": "\\frac{S}{k_B} = \\ln\\Omega",
                "Variables": ["kB", "Ω"],
                "Additional_Info": "Where kB is the Boltzman constant and Ω is the number of microstates."
            },

            # Partition Function
            {
                "Formula_Name": "Partition Function",
                "Unknown": "Z",
                "Formula": "Z = Σ_n e^(-βE_n)",
                "LaTeX": "Z = \\sum_n e^{-\\beta E_n}",
                "Variables": ["β", "En"],
                "Additional_Info": "Where β is equal to 1/kBT and En is the energy of state n."
            },
            {
                "Formula_Name": "Partition Function Using Density of States",
                "Unknown": "Z",
                "Formula": "Z = Σ_E g(E)e^(-βE)",
                "LaTeX": "Z = \\sum_E g(E)e^{-\\beta E}",
                "Variables": ["g(E)", "β", "E"],
                "Additional_Info": "Where g(E) is the density of states with energy E, β is equal to 1/kBT, and E is the energy level being summed over."
            },
            {
                "Formula_Name": "Partition Function (Continous Energy Spectrum)",
                "Unknown": "Z",
                "Formula": "Z = ∫g(E)e^(-βE)dE",
                "LaTeX": "Z = \\int g(E)e^{-\\beta E}dE",
                "Variables": ["g(E)", "β", "E", "dE"],
                "Additional_Info": "Where g(E) is the density of states with energy E, β is equal to 1/kBT, E is the energy level being summed over, and dE is the change in energy."
            },
            {
                "Formula_Name": "Partition Function",
                "Unknown": "Z",
                "Formula": "Z = e^(-βF)",
                "LaTeX": "Z = e^{-\\beta F}",
                "Variables": ["β", "F"],
                "Additional_Info": "Where β is equal to 1/kBT and F is the Helmholtz free energy."
            },
            {
                "Formula_Name": "RNA Partition Function",
                "Unknown": "Z",
                "Formula": "Z = Σ_s e^(-βE_s)",
                "LaTeX": "Z = \\sum_s e^{-\\beta E_s}",
                "Variables": ["β", "Es"],
                "Additional_Info": "Where β is equal to 1/kBT and Es is the energy of the specific microstate s."
            },

            # Energy Probability Distribution
            {
                "Formula_Name": "Energy Probability Distribution",
                "Unknown": "P(E)",
                "Formula": "P(E) = g(E)e^(-βE)/Z",
                "LaTeX": "P(E) = \\frac{g(E)e^{-\\beta E}}{Z}",
                "Variables": ["g(E)", "β", "E", "Z"],
                "Additional_Info": "Where g(E) is the density of states with energy E, β is equal to 1/kBT, E is the energy level being summed over, and Z is the partition function."
            },

            # Maxwell-Boltzmann Velocity Distribution
            {
                "Formula_Name": "Maxwell-Boltzmann Velocity Distribution",
                "Unknown": "P(vα)",
                "Formula": "P(v_α) = (m/2πk_BT)^(1/2)e^(-mv^2_α/2k_BT)",
                "LaTeX": "P(v_\\alpha) = \\left(\\frac{m}{2\\pi k_B T}\\right)^{1/2}e^{-\\frac{mv^{2}_{\\alpha}}{2k_B T}}",
                "Variables": ["m", "kB", "T", "vα"],
                "Additional_Info": "Where m is mass, kB is the Boltzman constant, T is temperature, and vα is the velocity component in the α direction (where α could represent x, y, or z coordinate)."
            },
            {
                "Formula_Name": "Maxwell-Boltzmann Velocity Distribution",
                "Unknown": "P(vα)",
                "Formula": "P(v_α) = sqrt(βm/2π)e^(-βmv^2_α/2)",
                "LaTeX": "P(v_\\alpha) = \\sqrt{\\frac{\\beta m}{2\\pi}}e^{-\\frac{\\beta mv^{2}_{\\alpha}}{2}}",
                "Variables": ["β", "m", "vα"],
                "Additional_Info": "Where β is equal to 1/kBT, m is the mass, and vα is the velocity component in the α direction (where α could represent x, y, or z coordinate)."
            },
            {
                "Formula_Name": "Maxwell-Boltzmann Speed Distribution",
                "Unknown": "P(v)",
                "Formula": "P(v) = 4π(m/2πk_BT)^(3/2)v^2e^(-mv^2/2k_BT)",
                "LaTeX": "P(v) = 4\\pi\\left(\\frac{m}{2\\pi k_B T}\\right)^{3/2}v^2e^{-\\frac{mv^2}{2k_B T}}",
                "Variables": ["m", "kB", "T", "v"],
                "Additional_Info": "Where m is mass, kB is the Boltzman constant, T is temperature, and v is velocity."
            },
            {
                "Formula_Name": "Most Probable Speed (Maxwell-Boltzmann Distribution)",
                "Unknown": "v_mp",
                "Formula": "v_mp = sqrt(2k_BT/m)",
                "LaTeX": "v_{mp} = \\sqrt{\\frac{2k_B T}{m}}",
                "Variables": ["kB", "T", "m"],
                "Additional_Info": "Where kB is the Boltzman constant, T is temperature, and m is particle mass."
            },

            # Average Kinetic Energy
            {
                "Formula_Name": "Average Kinetic Energy",
                "Unknown": "⟨E_kin⟩",
                "Formula": "⟨E_kin⟩ = 3k_BT/2",
                "LaTeX": "\\langle E_{kin}\\rangle = \\frac{3}{2}k_B T",
                "Variables": ["kB", "T"],
                "Additional_Info": "Where kB is the Boltzman constant and T is temperature."
            },
            {
                "Formula_Name": "Average Kinetic Energy",
                "Unknown": "⟨E_kin⟩",
                "Formula": "⟨E_kin⟩ = ⟨mv^2/2⟩",
                "LaTeX": "\\langle E_{kin}\\rangle = \\left\\langle\\frac{mv^2}{2}\\right\\rangle",
                "Variables": ["m", "v"],
                "Additional_Info": "Where m is mass and v is velocity."
            },
            {
                "Formula_Name": "Average Kinetic Energy",
                "Unknown": "⟨E_kin⟩",
                "Formula": "⟨E_kin⟩ = (m/2)(⟨v^2_x⟩+⟨v^2_y⟩+⟨v^2_z⟩)",
                "LaTeX": "\\langle E_{kin}\\rangle = \\frac{m}{2}(\\langle v^{2}_x\\rangle + \\langle v^{2}_y\\rangle + \\langle v^{2}_z\\rangle)",
                "Variables": ["m", "vx", "vy", "vz"],
                "Additional_Info": "Where m is mass, and vx, vy, and vz are velocities in different directions."
            },
            {
                "Formula_Name": "Average Kinetic Energy",
                "Unknown": "⟨E_kin⟩",
                "Formula": "⟨E_kin⟩ = 3m⟨v^2_α⟩/2",
                "LaTeX": "\\langle E_{kin}\\rangle = \\frac{3m}{2}\\langle v^{2}_{\\alpha}\\rangle",
                "Variables": ["m", "vα"],
                "Additional_Info": "Where m is mass and vα is the velocity component in the α direction (where α could represent x, y, or z coordinate)."
            },
            {
                "Formula_Name": "Mean Square Velocity Component",
                "Unknown": "⟨v^2_α⟩",
                "Formula": "⟨v^2_α⟩ = k_BT/m",
                "LaTeX": "\\langle v^{2}_{\\alpha}\\rangle = \\frac{k_B T}{m}",
                "Variables": ["kB", "T", "m"],
                "Additional_Info": "Where kB is the Boltzman constant, T is temperature, and m is particle mass."
            },

            # Rate of Change of Concentration
            {
                "Formula_Name": "Rate of Change of Concentration",
                "Unknown": "dx/dt",
                "Formula": "d[x]/dt = J_f - J_b",
                "LaTeX": "\\frac{d[x]}{dt} = J_f - J_b",
                "Variables": ["Jf", "Jb"],
                "Additional_Info": "Where Jf is the forward flux and Jb is the backward flux."
            },
            {
                "Formula_Name": "Rate of Change of Concentration",
                "Unknown": "dx/dt",
                "Formula": "d[x]/dt = k_f[A][B] - k_b[C]",
                "LaTeX": "\\frac{d[x]}{dt} = k_f[A][B] - k_b[C]",
                "Variables": ["kf", "kb", "A", "B", "C"],
                "Additional_Info": "Where kf is the forward rate constant, kb is the backward rate constant, and A, B, and C are concentrations."
            },
            {
                "Formula_Name": "Rate of Change of Concentration (For Reaction A + B ⇌ C)",
                "Unknown": "dA/dt",
                "Formula": "d[A]/dt = -k_f[A][B] + k_b[C]",
                "LaTeX": "\\frac{d[A]}{dt} = -k_f[A][B] + k_b[C]",
                "Variables": ["kf", "kb", "A", "B", "C"],
                "Additional_Info": "Where kf is the forward rate constant, kb is the backward rate constant, and A, B, and C are concentrations."
            },
            {
                "Formula_Name": "Rate of Change of Concentration (For Reaction A + B ⇌ C)",
                "Unknown": "dB/dt",
                "Formula": "d[B]/dt = -k_f[A][B] + k_b[C]",
                "LaTeX": "\\frac{d[B]}{dt} = -k_f[A][B] + k_b[C]",
                "Variables": ["kf", "kb", "A", "B", "C"],
                "Additional_Info": "Where kf is the forward rate constant, kb is the backward rate constant, and A, B, and C are concentrations."
            },
            {
                "Formula_Name": "Rate of Change of Concentration (For Reaction A + B ⇌ C)",
                "Unknown": "dC/dt",
                "Formula": "d[C]/dt = k_f[A][B] - k_b[C]",
                "LaTeX": "\\frac{d[C]}{dt} = k_f[A][B] - k_b[C]",
                "Variables": ["kf", "kb", "A", "B", "C"],
                "Additional_Info": "Where kf is the forward rate constant, kb is the backward rate constant, and A, B, and C are concentrations."
            },

            # Van't Hoff Equation
            {
                "Formula_Name": "Van't Hoff Equation",
                "Unknown": "Keq",
                "Formula": "K_eq = e^(-ΔG_0/RT)",
                "LaTeX": "K_{eq} = e^{-\\frac{\\Delta G_0}{RT}}",
                "Variables": ["ΔG0", "R", "T"],
                "Additional_Info": "Where ΔG0 is the standard Gibbs free energy change, R is the gas constant, and T is temperature."
            },
            {
                "Formula_Name": "Van't Hoff Equation",
                "Unknown": "Keq",
                "Formula": "K_eq = [C]^v_C[D]^v_D/([A]^v_A[B]^v_B)",
                "LaTeX": "K_{eq} = \\frac{[C]^{v_C}[D]^{v_D}}{[A]^{v_A}[B]^{v_B}}",
                "Variables": ["C", "D", "A", "B", "vC", "vD", "vA", "vB"],
                "Additional_Info": "Where A,B,C,D are concentrations of reactants and products, and vA,vB,vC,vD are stoichiometric coefficients."
            },
            {
                "Formula_Name": "Van't Hoff Equation",
                "Unknown": "Keq",
                "Formula": "K_eq = (J_b k_f)/(J_f k_b) = k_f/k_b",
                "LaTeX": "K_{eq} = \\frac{J_b k_f}{J_f k_b} = \\frac{k_f}{k_b}",
                "Variables": ["Jb", "kf", "Jf", "kb"],
                "Additional_Info": "Where Jf and Jb are forward and backward fluxes and kf and kb are forward and backward rate constants."
            },
            {
                "Formula_Name": "Van't Hoff Equation (Rearranged to Solve for ΔG_0)",
                "Unknown": "ΔG_0",
                "Formula": "ΔG_0 = -RT ln K_eq",
                "LaTeX": "\\Delta G_0 = -RT\\ln K_{eq}",
                "Variables": ["R", "T", "Keq"],
                "Additional_Info": "Where R is the gas constant, T is temperature. and Keq is the equilibrium constant."
            },
            {
                "Formula_Name": "Van't Hoff Equation (logarithmic form)",
                "Unknown": "ln(Keq)",
                "Formula": "ln K_eq = -ΔG_0/(RT)",
                "LaTeX": "\\ln K_{eq} = -\\frac{\\Delta G_0}{RT}",
                "Variables": ["ΔG0", "R", "T"],
                "Additional_Info": "Where ΔG0 is the standard Gibbs free energy change, R is the gas constant, and T is temperature."
            },
            {
                "Formula_Name": "Van't Hoff Equation (logarithmic form using ΔG₀ = ΔH₀ - TΔS₀)",
                "Unknown": "ln(Keq)",
                "Formula": "ln K_eq = -ΔH_0/(RT) + ΔS_0/R",
                "LaTeX": "\\ln K_{eq} = -\\frac{\\Delta H_0}{RT} + \\frac{\\Delta S_0}{R}",
                "Variables": ["ΔH0", "R", "T", "ΔS0"],
                "Additional_Info": "Where R is the gas constant, T is temperature, and ΔH0 and ΔS0 are changes in enthalpy and entropy."
            },

            # Arrhenius Equation
            {
                "Formula_Name": "Arrhenius Equation",
                "Unknown": "k",
                "Formula": "k = k_0 e^(-ΔE/(k_B T))",
                "LaTeX": "k = k_0 e^{-\\frac{\\Delta E}{k_B T}}",
                "Variables": ["k0", "ΔE", "kB", "T"],
                "Additional_Info": "Where k0 is the pre-exponential factor, ΔE is the activation energy, Kb is the Boltzmann constant, and T is temperature."
            },
            {
                "Formula_Name": "Arrhenius Equation",
                "Unknown": "k",
                "Formula": "k = k_0 e^(-ΔG^‡/(RT))",
                "LaTeX": "k = k_0 e^{-\\frac{\\Delta G^{\\ddagger}}{RT}}",
                "Variables": ["k0", "ΔG‡", "R", "T"],
                "Additional_Info": "Where k0 is the pre-exponential factor, ΔG‡ is the Gibbs free energy of activation, R is the gas constant, and T is temperature."
            },
            {
                "Formula_Name": "Arrhenius Equation (Using the Gas Constant R and Activation Energy Ea)",
                "Unknown": "k",
                "Formula": "k = A e^(-E_a/(RT))",
                "LaTeX": "k = A e^{-\\frac{E_a}{RT}}",
                "Variables": ["A", "Ea", "R", "T"],
                "Additional_Info": "Where A is the frequency or pre-exponential factor and e^(-Ea/RT) represents the fraction of collisions that have enough energy to overcome the activation barrier."
            },
            {
                "Formula_Name": "Arrhenius Equation (Using Activation Energy Ea)",
                "Unknown": "k",
                "Formula": "k = k_0 e^(-E_a/(k_B T))",
                "LaTeX": "k = k_0 e^{-\\frac{E_a}{k_B T}}",
                "Variables": ["k0", "Ea", "kB", "T"],
                "Additional_Info": "Where k0 is the pre-exponential factor, Ea is the energy of activation, kB is the Boltzman constant, and T is temperature."
            },
            {
                "Formula_Name": "Arrhenius Equation (Logarithmic Form)",
                "Unknown": "ln(k)",
                "Formula": "ln(k) = ln(k_0) - ΔE/(k_B T)",
                "LaTeX": "\\ln(k) = \\ln(k_0) - \\frac{\\Delta E}{k_B T}",
                "Variables": ["k0", "ΔE", "kB", "T"],
                "Additional_Info": "Where k0 is the pre-exponential factor, ΔE is the activation energy, Kb is the Boltzmann constant, and T is temperature."
            },
            {
                "Formula_Name": "Arrhenius Equation (Ratio Form)",
                "Unknown": "k1/k2",
                "Formula": "k_1/k_2 = e^(-(ΔE/k_B)(1/T_2 - 1/T_1))",
                "LaTeX": "\\frac{k_1}{k_2} = e^{-\\frac{\\Delta E}{k_B}\\left(\\frac{1}{T_2}-\\frac{1}{T_1}\\right)}",
                "Variables": ["ΔE", "kB", "T1", "T2"],
                "Additional_Info": "Where ΔE is the activation energy, kB is the Boltzman constant, and T1 and T2 are temperatures."
            },

            # Henderson-Hasselbalch Equation
            {
                "Formula_Name": "Henderson-Hasselbalch Equation (Solving for pK)",
                "Unknown": "pK",
                "Formula": "pK = pH - log([A^-]/[HA])",
                "LaTeX": "pK = pH - \\log\\left(\\frac{[A^-]}{[HA]}\\right)",
                "Variables": ["pH", "A-", "HA"],
                "Additional_Info": "Where pH is the negative log of hydrogen ion concentration, and ([A⁻]/[HA]) is the ratio of deprotonated to protonated acid."
            },
            {
                "Formula_Name": "Henderson-Hasselbalch Equation (Solving for pH)",
                "Unknown": "pH",
                "Formula": "pH = pKa + log([A^-]/[HA])",
                "LaTeX": "pH = pKa + \\log\\left(\\frac{[A^-]}{[HA]}\\right)",
                "Variables": ["pKa", "A-", "HA"],
                "Additional_Info": "Where pKa is the negative log of the equilibrium constant and ([A⁻]/[HA]) is the ratio of deprotonated to protonated acid."
            },
            {
                "Formula_Name": "Henderson-Hasselbalch Equation",
                "Unknown": "log(A-/HA)",
                "Formula": "log([A^-]/[HA]) = pH - pKa",
                "LaTeX": "\\log\\left(\\frac{[A^-]}{[HA]}\\right) = pH - pKa",
                "Variables": ["pH", "pKa"],
                "Additional_Info": "Where pH is the negative log of hydrogen ion concentration and pKa is the negative log of the equilibrium constant."
            },
            {
                "Formula_Name": "Henderson-Hasselbalch Equation",
                "Unknown": "A-/HA",
                "Formula": "[A^-]/[HA] = 10^(pH-pKa)",
                "LaTeX": "\\frac{[A^-]}{[HA]} = 10^{pH-pKa}",
                "Variables": ["pH", "pKa"],
                "Additional_Info": "Where pH is the negative log of hydrogen ion concentration and pKa is the negative log of the equilibrium constant."
            },

            # Mean Square Displacement in Diffusion
            {
                "Formula_Name": "Mean Square Displacement in Diffusion (1D)",
                "Unknown": "⟨r^2⟩",
                "Formula": "⟨r^2⟩ = 2Dt (in 1D)",
                "LaTeX": "\\langle r^2 \\rangle = 2Dt \\text{ (in 1D)}",
                "Variables": ["D", "t"],
                "Additional_Info": "Where D is the diffusion coefficient and t is time."
            },
            {
                "Formula_Name": "Mean Square Displacement in Diffusion (2D)",
                "Unknown": "⟨r^2⟩",
                "Formula": "⟨r^2⟩ = 4Dt (in 2D)",
                "LaTeX": "\\langle r^2 \\rangle = 4Dt \\text{ (in 2D)}",
                "Variables": ["D", "t"],
                "Additional_Info": "Where D is the diffusion coefficient and t is time."
            },
            {
                "Formula_Name": "Mean Square Displacement in Diffusion (3D)",
                "Unknown": "⟨r^2⟩",
                "Formula": "⟨r^2⟩ = 6Dt",
                "LaTeX": "\\langle r^2 \\rangle = 6Dt",
                "Variables": ["D", "t"],
                "Additional_Info": "Where D is the diffusion coefficient and t is time."
            },
            {
                "Formula_Name": "Root Mean Square Displacement in Diffusion (3D)",
                "Unknown": "sqrt{⟨r^2⟩}",
                "Formula": "sqrt{⟨r^2⟩} = sqrt{6Dt}",
                "LaTeX": "\\sqrt{\\langle r^2 \\rangle} = \\sqrt{6Dt}",
                "Variables": ["D", "t"],
                "Additional_Info": "Where D is the diffusion coefficient and t is time."
            },

            # Fick's First Law of Diffusion
            {
                "Formula_Name": "Fick's First Law of Diffusion",
                "Unknown": "J",
                "Formula": "J = -D(∂C/∂x)",
                "LaTeX": "J = -D\\frac{\\partial C}{\\partial x}",
                "Variables": ["D", "∂C/∂x"],
                "Additional_Info": "Where D is the diffusion coefficient and (∂C/∂x) is a concentration gradient."
            },
            {
                "Formula_Name": "Fick's First Law of Diffusion",
                "Unknown": "J",
                "Formula": "J = (1/2)fb(C_1-C_2)",
                "LaTeX": "J = \\frac{1}{2}fb(C_1-C_2)",
                "Variables": ["fb", "C1", "C2"],
                "Additional_Info": ""
            },

            # Langevin Equation
            {
                "Formula_Name": "Langevin Equation",
                "Unknown": "m(d^2r/dt^2)",
                "Formula": "m(d^2r/dt^2) = F - ξ(dr/dt) + R",
                "LaTeX": "m\\frac{d^2r}{dt^2} = F - \\xi\\frac{dr}{dt} + R",
                "Variables": ["F", "ξ", "dr/dt", "R"],
                "Additional_Info": "Where F is external force, ξ is the friction coefficient, dr/dt is velocity, and R is random force."
            },
            {
                "Formula_Name": "Langevin Equation (Using Velocity)",
                "Unknown": "ma",
                "Formula": "ma = F - ξv + R",
                "LaTeX": "ma = F - \\xi v + R",
                "Variables": ["F", "ξ", "v", "R"],
                "Additional_Info": "Where F is external force, v is velocity, R is random force, and ξ is the friction coefficient."
            },

            # Stokes Law
            {
                "Formula_Name": "Stokes Law",
                "Unknown": "ξ",
                "Formula": "ξ = 6πηa",
                "LaTeX": "\\xi = 6\\pi\\eta a",
                "Variables": ["η", "a"],
                "Additional_Info": "Where η is fluid viscosity and a is particle radius."
            },
            {
                "Formula_Name": "Stokes Law (For a Sphere at a Liquid-Gas Interface)",
                "Unknown": "ξ",
                "Formula": "ξ = 4πηa",
                "LaTeX": "\\xi = 4\\pi\\eta a",
                "Variables": ["η", "a"],
                "Additional_Info": "Where η is fluid viscosity and a is particle radius."
            },
            {
                "Formula_Name": "Damping Coefficient per Unit Mass",
                "Unknown": "γ",
                "Formula": "γ = 6πηa/m",
                "LaTeX": "\\gamma = \\frac{6\\pi\\eta a}{m}",
                "Variables": ["η", "a", "m"],
                "Additional_Info": "Where η is fluid viscosity, a is particle radius, and m is particle mass."
            },

            # Reynolds Number
            {
                "Formula_Name": "Reynolds Number",
                "Unknown": "R",
                "Formula": "R = vr/ν",
                "LaTeX": "R = \\frac{vr}{\\nu}",
                "Variables": ["v", "r", "ν"],
                "Additional_Info": "Where v is velocity, r is characteristic length/radius, and ν is kinematic viscosity."
            },
            {
                "Formula_Name": "Reynolds Number (using friction coefficient)",
                "Unknown": "R",
                "Formula": "R = mv/(ξr)",
                "LaTeX": "R = \\frac{mv}{\\xi r}",
                "Variables": ["m", "v", "ξ", "r"],
                "Additional_Info": "Where m is mass, v is velocity, ξ is friction coefficient, and r is radius."
            },
            {
                "Formula_Name": "Reynolds Number (using dynamic viscosity and density)",
                "Unknown": "R",
                "Formula": "R = vr/(η/ρ)",
                "LaTeX": "R = \\frac{vr}{\\eta/\\rho}",
                "Variables": ["v", "r", "η", "ρ"],
                "Additional_Info": "Where v is velocity, r is radius, η is viscosity, and ρ is density."
            },

            # Schrödinger Equation
            {
                "Formula_Name": "Time-Independent Schrödinger Equation",
                "Unknown": "Hψ",
                "Formula": "Ĥψ = Eψ",
                "LaTeX": "\\hat{H}\\psi = E\\psi",
                "Variables": ["E", "ψ"],
                "Additional_Info": "Where Ĥ is the Hamiltonian operator representing the total energy operator, E is energy eigenvalues, and ψ is the wave function."
            },
            {
                "Formula_Name": "Time-Independent Schrödinger Equation",
                "Unknown": "Eψ",
                "Formula": "(-h^2/(8π^2m)∇^2 + E_pot)ψ = Eψ",
                "LaTeX": "\\left(-\\frac{h^2}{8\\pi^2m}\\nabla^2 + E_{pot}\\right)\\psi = E\\psi",
                "Variables": ["h", "m", "∇^2", "Epot", "ψ", "E"],
                "Additional_Info": "Where h is Planck's constant, m is mass, ∇^2 is the Laplacian operator, Epot is potential energy, ψ is the wave function, and E is energy eigenvalues."
            },

            # Hydrogen Atom Energy Levels
            {
                "Formula_Name": "Hydrogen Atom Energy Levels",
                "Unknown": "En",
                "Formula": "E_n = -2π^2mq^4/(h^2n^2)",
                "LaTeX": "E_n = -\\frac{2\\pi^2mq^4}{h^2n^2}",
                "Variables": ["m", "q", "h", "n"],
                "Additional_Info": "Where m is electron mass, q is elementary charge, h is Planck's constant, and n is principal quantum number."
            },

            # Ionic Bond Energy
            {
                "Formula_Name": "Ionic Bond Energy",
                "Unknown": "E_IB",
                "Formula": "E_IB = q^2/(4πε_0r_0)",
                "LaTeX": "E_{IB} = \\frac{q^2}{4\\pi\\varepsilon_0 r_0}",
                "Variables": ["q", "ε0", "r0"],
                "Additional_Info": "Where q is the charge of ions, ε0 is the permittivity of free space, and r0 is the equilibrium distance between ions."
            },
            {
                "Formula_Name": "Ionic Bond Energy",
                "Unknown": "E_IB",
                "Formula": "E_IB = W = ∫Fdr",
                "LaTeX": "E_{IB} = W = \\int F\\,dr",
                "Variables": ["F", "dr"],
                "Additional_Info": "Where F is force and dr is the differential displacement."
            },
            {
                "Formula_Name": "Ionic Bond Energy (Coulomb's Law)",
                "Unknown": "F",
                "Formula": "F = q^2/(4πε_0r^2)",
                "LaTeX": "F = \\frac{q^2}{4\\pi\\varepsilon_0 r^2}",
                "Variables": ["q", "ε0", "r"],
                "Additional_Info": "Where q is the charge of ions, ε0 is the permittivity of free space, and r is the distance between ions."
            },

            {"Formula_Name": "Lennard-Jones Potential", "Unknown": "V_LJ",

             "Formula": "V_LJ(r) = ε_h[((σ/r)^12)-2((σ/r)^6)]",
               "LaTeX": "V_{LJ}(r) = \\varepsilon_h[(\\frac{\\sigma}{r})^{12}-2(\\frac{\\sigma}{r})^6]",
                "Variables": ["ε_h", "σ", "r"], "Additional_Info": "Where ε_h is the well depth (energy minimum), σ is the collision diameter (distance at which potential is zero), and r is the interatomic distance."
            },
            {
  "Formula_Name": "Lennard-Jones Potential (Alternative Form)",
  "Unknown": "V_LJ",
  "Formula": "V_LJ(r) = 4ε[((σ/r)^12)-((σ/r)^6)]",
  "LaTeX": "V_{LJ}(r) = 4\\varepsilon[(\\frac{\\sigma}{r})^{12}-(\\frac{\\sigma}{r})^6]",
  "Variables": ["ε", "σ", "r"],
  "Additional_Info": "Where ε is the well depth, σ is the collision diameter, and r is the interatomic distance."
},
            {
  "Formula_Name": "Lennard-Jones Potential (Using Position of Minimum)",
  "Unknown": "V_LJ",
  "Formula": "V_LJ(r) = ε[((r_m/r)^12)-2((r_m/r)^6)]",
  "LaTeX": "V_{LJ}(r) = \\varepsilon[(\\frac{r_m}{r})^{12}-2(\\frac{r_m}{r})^6]",
  "Variables": ["ε", "r_m", "r"],
  "Additional_Info": "Where ε is the well depth, r_m is the position of minimum, and r is the interatomic distance."
},
            {
  "Formula_Name": "Force from Lennard-Jones Potential",
  "Unknown": "f",
  "Formula": "f(r) = (12ε_h/r)[((σ/r)^12)-((σ/r)^6)]",
  "LaTeX": "f(r) = \\frac{12\\varepsilon_h}{r}[(\\frac{\\sigma}{r})^{12}-(\\frac{\\sigma}{r})^6]",
  "Variables": ["ε_h", "σ", "r"],
  "Additional_Info": "Where ε_h is the well depth (energy minimum), σ is the collision diameter, and r is the interatomic distance."
},
            {
  "Formula_Name": "Force from Lennard-Jones Potential (Derivative Form)",
  "Unknown": "f",
  "Formula": "f(r) = -∂V/∂r",
  "LaTeX": "f(r) = -\\frac{\\partial V}{\\partial r}",
  "Variables": ["V", "r"],
  "Additional_Info": "Where ∂V/∂r is the partial derivative of potential energy with respect to interatomic distance."
},
            {
  "Formula_Name": "Electrostatic Energy (Coulomb's Potential)",
  "Unknown": "E_EL",
  "Formula": "E_EL = (q_1q_2)/(4πε_0r)",
  "LaTeX": "E_{EL} = \\frac{q_1q_2}{4\\pi\\varepsilon_0r}",
  "Variables": ["q1", "q2", "ε0", "r"],
  "Additional_Info": "Where q1 and q2 are the charges of two particles, ε0 is the permittivity of free space, and r is the distance between the charges."
},
            {
  "Formula_Name": "Bond-Length Potential",
  "Unknown": "V_bl",
  "Formula": "V_bl(r) = k_α(r-r_0)^2",
  "LaTeX": "V_{bl}(r) = k_{\\alpha}(r-r_0)^2",
  "Variables": ["kα", "r", "r0"],
  "Additional_Info": "Where kα is the bond force constant, r is the current bond length, and r0 is the equilibrium bond length."
},
            {
  "Formula_Name": "Bond-Angle Potential",
  "Unknown": "V_ba",
  "Formula": "V_ba(θ) = k_θ(θ-θ_0)^2",
  "LaTeX": "V_{ba}(\\theta) = k_{\\theta}(\\theta-\\theta_0)^2",
  "Variables": ["kθ", "θ", "θ0"],
  "Additional_Info": "Where kθ is the angle force constant, θ is the current bond angle, and θ0 is the equilibrium bond angle."
},
            {
  "Formula_Name": "Dihedral Angle Potential",
  "Unknown": "V",
  "Formula": "V(ϕ) = (V_0/2)(1+cos(ϕ_0+nϕ))",
  "LaTeX": "V(\\phi) = \\frac{V_0}{2}(1+\\cos(\\phi_0+n\\phi))",
  "Variables": ["V0", "ϕ0", "n", "ϕ"],
  "Additional_Info": "Where V0 is the barrier height, ϕ0 is the phase shift, n is the periodicity, and ϕ is the dihedral angle."
},
            {
  "Formula_Name": "Improper Dihedral Angle Potential",
  "Unknown": "V_ia",
  "Formula": "V_ia(ψ) = k_ψ(ψ-ψ_0)^2",
  "LaTeX": "V_{ia}(\\psi) = k_{\\psi}(\\psi-\\psi_0)^2",
  "Variables": ["kψ", "ψ", "ψ0"],
  "Additional_Info": "Where kψ is the improper dihedral force constant, ψ is the current improper dihedral angle, and ψ0 is the equilibrium improper dihedral angle."
},
            {
  "Formula_Name": "Total Molecular Energy Function",
  "Unknown": "E_tot",
  "Formula": "E_tot = Σ_{all BL}k_a(r-r_0)^2 + Σ_{all BA}k_θ(θ-θ_0)^2 + Σ_{all DA}(V_0/2)(1+cos(ϕ_0+nϕ)) + Σ_{all IA}k_ψ(ψ-ψ_0)^2 + Σ_{ij}(q_iq_j)/(4πε_0r) + Σ_{ij}ε_h[((σ/r_{ij})^12)-2((σ/r_{ij})^6)]",
  "LaTeX": "E_{tot} = \\Sigma_{allBL}k_a(r-r_0)^2+\\Sigma_{allBA}k_{\\theta}(\\theta-\\theta_0)^2+\\Sigma_{allDA}\\frac{V_0}{2}(1+\\cos(\\phi_0+n\\phi))+\\Sigma_{allIA}k_{\\psi}(\\psi-\\psi_0)^2+\\Sigma_{ij}\\frac{q_iq_j}{4\\pi\\varepsilon_0r}+\\Sigma_{ij}\\varepsilon_h[(\\frac{\\sigma}{r_{ij}})^{12}-2(\\frac{\\sigma}{r_{ij}})^6]",
  "Variables": ["ka", "r", "r0", "kθ", "θ", "θ0", "V0", "ϕ0", "n", "ϕ", "kψ", "ψ", "ψ0", "qi", "qj", "ε0", "εh", "σ", "rij"],
  "Additional_Info": "Where all terms include force constants, geometry parameters, equilibrium values, charges, and Lennard-Jones parameters, summing bond lengths, bond angles, dihedral angles, improper dihedrals, electrostatic interactions, and van der Waals forces."
},
            {
  "Formula_Name": "Salt Bridge Energy",
  "Unknown": "E_SB",
  "Formula": "E_SB = (1/(4πε_0))(q_1q_2/r_0)",
  "LaTeX": "E_{SB} = \\frac{1}{4\\pi\\varepsilon_0}\\frac{q_1q_2}{r_0}",
  "Variables": ["q1", "q2", "r0", "ε0"],
  "Additional_Info": "Where q1 and q2 are the charges of interacting groups, r0 is the distance between charges, and ε0 is the permittivity of free space."
},
            {
  "Formula_Name": "Salt Bridge Energy (using Coulomb's constant)",
  "Unknown": "E_SB",
  "Formula": "E_SB = (k_eq_1q_2)/r_0",
  "LaTeX": "E_{SB} = \\frac{k_eq_1q_2}{r_0}",
  "Variables": ["ke", "q1", "q2", "r0"],
  "Additional_Info": "Where ke is Coulomb's constant, q1 and q2 are the charges of interacting groups, and r0 is the distance between charges."
},
            {
  "Formula_Name": "ATP Energy",
  "Unknown": "E_ATP",
  "Formula": "E_ATP = (1/(4πε_0ε))(q_1q_2/r_0)",
  "LaTeX": "E_{ATP} = \\frac{1}{4\\pi\\varepsilon_0\\varepsilon}\\frac{q_1q_2}{r_0}",
  "Variables": ["ε", "ε0", "q1", "q2", "r0"],
  "Additional_Info": "Where ε0 is the permittivity of free space, ε is the dielectric constant of the medium, q1 and q2 are the charges, and r0 is the distance between charges."
},
            {
  "Formula_Name": "Dielectric Constant formula",
  "Unknown": "ε",
  "Formula": "ε = (1/(4πε_0))(q_1q_2/(r_0E))",
  "LaTeX": "\\varepsilon = \\frac{1}{4\\pi\\varepsilon_0}\\frac{q_1q_2}{r_0E}",
  "Variables": ["ε0", "q1", "q2", "r0", "E"],
  "Additional_Info": "Where ε0 is the permittivity of free space, q1 and q2 are the charges, r0 is the distance between charges, and E is the energy."
},
            {
  "Formula_Name": "Radial Distribution Function",
  "Unknown": "g(r)",
  "Formula": "g(r) = ⟨h(r)⟩/V(r)",
  "LaTeX": "g(r) = \\frac{\\langle h(r)\\rangle}{V(r)}",
  "Variables": ["⟨h(r)⟩", "V(r)"],
  "Additional_Info": "Where ⟨h(r)⟩ is the average histogram of particles at distance r, and V(r) is the volume of the concentric layer at distance r."
},
            {
  "Formula_Name": "Ionization Energy",
  "Unknown": "E_I",
  "Formula": "E_I = hv = hc/λ",
  "LaTeX": "E_I = h\\nu = \\frac{hc}{\\lambda}",
  "Variables": ["h", "v", "λ", "c"],
  "Additional_Info": "Where h is Planck's constant, v is the frequency of radiation, λ is the wavelength of radiation, and c is the speed of light."
},
            {
  "Formula_Name": "Surface-to-Volume Ratio Formula",
  "Unknown": "Ns/N",
  "Formula": "Ns/N = (ρ_wV_S)/(ρ_wV) = (4πr^2d)/((4/3)πr^3) ∼ d/r",
  "LaTeX": "\\frac{N_S}{N}=\\frac{\\rho_wV_S}{\\rho_wV} = \\frac{4\\pi r^2d}{\\frac{4}{3}\\pi r^3}\\sim\\frac{d}{r}",
  "Variables": ["r", "d"],
  "Additional_Info": "Where r is the radius of the sphere, and d is the thickness of the surface layer."
},
            {
  "Formula_Name": "Surface-to-Volume Ratio (simplified form)",
  "Unknown": "Ns/N",
  "Formula": "Ns/N = 3d/r",
  "LaTeX": "\\frac{N_S}{N} = \\frac{3d}{r}",
  "Variables": ["d", "r"],
  "Additional_Info": "Where d is the thickness of the surface layer, and r is the radius of the sphere."
},
            {
  "Formula_Name": "Boltzmann Distribution of Ion Concentration",
  "Unknown": "Ci",
  "Formula": "C_i = C_{i,0}e^(-z_iqeφ(r)/(k_BT))",
  "LaTeX": "C_i = C_{i,0}e^{-\\frac{z_iq_e\\phi(r)}{k_BT}}",
  "Variables": ["Ci0", "zi", "qe", "φ(r)", "kB", "T"],
  "Additional_Info": "Where Ci,0 is the bulk ion concentration, zi is the ion valence, qe is the elementary charge, φ(r) is the electric potential, kB is the Boltzmann constant, and T is the temperature."
},
            {
  "Formula_Name": "Boltzmann Distribution (energy form)",
  "Unknown": "Ci",
  "Formula": "C_i = C_{i,0}e^(-E_i/(k_BT))",
  "LaTeX": "C_i = C_{i,0}e^{-\\frac{E_i}{k_BT}}",
  "Variables": ["Ci0", "Ei", "kB", "T"],
  "Additional_Info": "Where Ci,0 is the bulk ion concentration, Ei is the energy of the ion, kB is the Boltzmann constant, and T is the temperature."
},
            {
  "Formula_Name": "Poisson Equation for Electrostatics",
  "Unknown": "∇^2φ",
  "Formula": "∇^2φ = -(1/(εε_0))ρ(r)",
  "LaTeX": "\\nabla^2\\phi = -\\frac{1}{\\varepsilon\\varepsilon_0}\\rho(r)",
  "Variables": ["ρ(r)", "ε0", "ε"],
  "Additional_Info": "Where ρ(r) is the charge density, ε0 is the permittivity of free space, and ε is the dielectric constant."
},
            {
  "Formula_Name": "Debye-Hückel Length",
  "Unknown": "rD",
  "Formula": "r_D = sqrt((εε_0k_BT)/(2FC_0q_e))",
  "LaTeX": "r_D = \\sqrt{\\frac{\\varepsilon\\varepsilon_0k_BT}{2FC_0q_e}}",
  "Variables": ["ε", "ε0", "kB", "T", "F", "C0", "qe"],
  "Additional_Info": "Where ε is the dielectric constant, ε0 is the permittivity of free space, kB is the Boltzmann constant, T is the temperature, F is the Faraday constant, C0 is the bulk salt concentration, and qe is the elementary charge."
},
            {
  "Formula_Name": "Einstein Relation/Diffusion-Mobility Relationship",
  "Unknown": "D",
  "Formula": "D = μ(k_BT/q)",
  "LaTeX": "D = \\mu\\frac{k_BT}{q}",
  "Variables": ["μ", "kB", "T", "q"],
  "Additional_Info": "Where μ is the mobility, kB is the Boltzmann constant, T is the temperature, and q is the charge."
},
            {
  "Formula_Name": "Electric Field from Potential Gradient",
  "Unknown": "E",
  "Formula": "E = -dφ/dx",
  "LaTeX": "E = -\\frac{d\\phi}{dx}",
  "Variables": ["dφ/dx"],
  "Additional_Info": "Where dφ/dx is the electric potential gradient."
},
            {
  "Formula_Name": "Electric Field from Potential Gradient (Vector Form)",
  "Unknown": "E",
  "Formula": "E = -∇φ",
  "LaTeX": "E = -\\nabla\\phi",
  "Variables": ["∇φ"],
  "Additional_Info": "Where ∇φ is the gradient of the electric potential."
},
            {
  "Formula_Name": "Nernst Equation",
  "Unknown": "ΔV",
  "Formula": "ΔV = -(k_BT/q)ln(C_in/C_out)",
  "LaTeX": "\\Delta V = -\\frac{k_BT}{q}\\ln\\frac{C_{in}}{C_{out}}",
  "Variables": ["kB", "T", "q", "Cin", "Cout"],
  "Additional_Info": "Where kB is the Boltzmann constant, T is the temperature, q is the charge, Cin is the ion concentration inside, and Cout is the ion concentration outside."
},
            {
  "Formula_Name": "Weighted Average Membrane Potential",
  "Unknown": "V_tot",
  "Formula": "V_tot = (Σ_ig_iΔV_i)/(Σ_ig_i)",
  "LaTeX": "V_{tot} = \\frac{\\sum_i g_i\\Delta V_i}{\\sum_i g_i}",
  "Variables": ["gi", "ΔVi"],
  "Additional_Info": "Where gi are the individual conductances, and ΔVi are the individual equilibrium potentials."
},
            {
  "Formula_Name": "Goldman-Hodgkin-Katz (GHK) Voltage Equation",
  "Unknown": "Vm",
  "Formula": "V_m = (RT/F)ln((g_K[K^+]_o + g_Na[Na^+]_o + g_Cl[Cl^-]_i)/(g_K[K^+]_i + g_Na[Na^+]_i + g_Cl[Cl^-]_o))",
  "LaTeX": "V_m = \\frac{RT}{F}\\ln(\\frac{g_K[K^+]_o + g_{Na}[Na^+]_o + g_{Cl}[Cl^-]_i}{g_K[K^+]_i + g_{Na}[Na^+]_i + g_{Cl}[Cl^-]_o})",
  "Variables": ["R", "T", "F", "gK", "gNa", "gCl", "[K+]i", "[Na+]i", "[Cl−]i", "[K+]o", "[Na+]o", "[Cl−]o"],
  "Additional_Info": "Where R is the gas constant, T is the temperature, F is the Faraday constant, gK, gNa, gCl are the ion conductances, and [K+]i, [Na+]i, [Cl−]i, [K+]o, [Na+]o, [Cl−]o are the ion concentrations inside and outside."
},
            {
  "Formula_Name": "Goldman-Hodgkin-Katz (GHK) Voltage Equation (Permeability Form)",
  "Unknown": "Vm",
  "Formula": "V_m = (RT/F)ln((P_K[K^+]_o + P_Na[Na^+]_o + P_Cl[Cl^-]_i)/(P_K[K^+]_i + P_Na[Na^+]_i + P_Cl[Cl^-]_o))",
  "LaTeX": "V_m = \\frac{RT}{F}\\ln(\\frac{P_K[K^+]_o + P_{Na}[Na^+]_o + P_{Cl}[Cl^-]_i}{P_K[K^+]_i + P_{Na}[Na^+]_i + P_{Cl}[Cl^-]_o})",
  "Variables": ["R", "T", "F", "PK", "PNa", "PCl", "[K+]i", "[Na+]i", "[Cl−]i", "[K+]o", "[Na+]o", "[Cl−]o"],
  "Additional_Info": "Where R is the gas constant, T is the temperature, F is the Faraday constant, PK, PNa, PCl are the ion permeabilities, and [K+]i, [Na+]i, [Cl−]i, [K+]o, [Na+]o, [Cl−]o are the ion concentrations inside and outside."
},
            {
  "Formula_Name": "Ohm's Law for Axial Current",
  "Unknown": "-I",
  "Formula": "-I = ΔV/R^O_i = (1/R_i)(dV/dx)",
  "LaTeX": "-I = \\frac{\\Delta V}{R^{O}_{i}} = \\frac{1}{R_i} \\frac{dV}{dx}",
  "Variables": ["ΔV", "RiO", "Ri", "dV/dx"],
  "Additional_Info": "Where ΔV is the voltage difference, RiO is the segment resistance, Ri is the resistance per unit length, and dV/dx is the voltage gradient."
},
            {
  "Formula_Name": "Ohm's Law for Axial Current (Discrete Form)",
  "Unknown": "I",
  "Formula": "I = -ΔV/(R_iΔx)",
  "LaTeX": "I = -\\frac{\\Delta V}{R_i\\Delta x}",
  "Variables": ["ΔV", "Ri", "Δx"],
  "Additional_Info": "Where ΔV is the voltage difference, Ri is the resistance per unit length, and Δx is the distance increment."
},
            {
  "Formula_Name": "Current Loss Equation",
  "Unknown": "−ΔI",
  "Formula": "−ΔI = V/R^O_m",
  "LaTeX": "−\\Delta I = \\frac{V}{R^{O}_{m}}",
  "Variables": ["V", "RmO"],
  "Additional_Info": "Where V is the membrane potential, and RmO is the membrane resistance."
},
            {
  "Formula_Name": "Current Loss Equation (Differential Form)",
  "Unknown": "-dI/dx",
  "Formula": "-dI/dx = V/R_m",
  "LaTeX": "-\\frac{dI}{dx} = \\frac{V}{R_m}",
  "Variables": ["V", "Rm"],
  "Additional_Info": "Where V is the membrane potential, and Rm is the membrane resistance per unit length."
},
            {
  "Formula_Name": "Current Loss Equation (Discrete Form)",
  "Unknown": "ΔI",
  "Formula": "ΔI = -VΔx/R_m",
  "LaTeX": "\\Delta I = -\\frac{V\\Delta x}{R_m}",
  "Variables": ["V", "Δx", "Rm"],
  "Additional_Info": "Where V is the membrane potential, Δx is the distance increment, and Rm is the membrane resistance per unit length."
},
            {
  "Formula_Name": "Cable Equation",
  "Unknown": "V",
  "Formula": "(R_m/R_i)(d^2V/dx^2) = V",
  "LaTeX": "\\frac{R_m}{R_i}\\frac{d^2V}{dx^2} = V",
  "Variables": ["Rm", "Ri", "x"],
  "Additional_Info": "Where Rm is the membrane resistance per unit length, Ri is the intracellular resistance per unit length, and x is the position along the axon."
},
            {
  "Formula_Name": "Cable Equation (With Space Constant)",
  "Unknown": "V",
  "Formula": "λ^2(d^2V/dx^2) = V",
  "LaTeX": "\\lambda^2\\frac{d^2V}{dx^2} = V",
  "Variables": ["λ", "x"],
  "Additional_Info": "Where λ is the space constant defined as sqrt(Rm/Ri), and x is the position along the axon."
},
            {
  "Formula_Name": "Cable Equation (Time Dependent)",
  "Unknown": "V",
  "Formula": "τ(∂V/∂t) = λ^2(d^2V/dx^2) - V",
  "LaTeX": "\\tau\\frac{\\partial V}{\\partial t} = \\lambda^2\\frac{d^2V}{dx^2} - V",
  "Variables": ["τ", "λ", "t", "x"],
  "Additional_Info": "Where τ is the membrane time constant, λ is the space constant, t is time, and x is position."
},
            {
  "Formula_Name": "Cable Equation (Steady State)",
  "Unknown": "V(x)",
  "Formula": "V(x) = V_0exp[-x/λ]",
  "LaTeX": "V(x) = V_0\\exp[-\\frac{x}{\\lambda}]",
  "Variables": ["V0", "x", "λ"],
  "Additional_Info": "Where V0 is the initial voltage, x is the distance from the voltage source, and λ is the space constant."
},
            {
  "Formula_Name": "Activation Variable Differential Equation",
  "Unknown": "dm/dt",
  "Formula": "dm/dt = α_m(1-m) - β_mm",
  "LaTeX": "\\frac{dm}{dt} = \\alpha_m(1-m) - \\beta_m m",
  "Variables": ["αm", "βm", "m"],
  "Additional_Info": "Where αm is the activation rate constant, βm is the deactivation rate constant, and m is the current activation state."
},
            {
  "Formula_Name": "Inactivation Variable Differential Equation",
  "Unknown": "dh/dt",
  "Formula": "dh/dt = α_h(1-h) - β_hh",
  "LaTeX": "\\frac{dh}{dt} = \\alpha_h(1-h) - \\beta_h h",
  "Variables": ["αh", "βh", "h"],
  "Additional_Info": "Where αh is the inactivation rate constant, βh is the removal of inactivation rate constant, and h is the current inactivation state."
},
    {
"Formula_Name": "Frequency of Neutron",
"Unknown": "v",
"Formula": "v = \frac{λ}{τ} = \frac{1}{C_m\sqrt{R_m}}",
"LaTeX": "v = \\frac{\\lambda}{\\tau} = \\frac{1}{C_m\\sqrt{R_m}}",
"Variables": ["λ", "τ", "Cm", "Rm"],
"Additional_Info": "Where v is frequency, λ is wavelength, τ is period, Cm is capacitance, and Rm is resistance."
},
    {
"Formula_Name": "Chemical Potential of Component i",
"Unknown": "μi",
"Formula": "μ_i = μ_i^0 + RT ln(f_iX_i)",
"LaTeX": "\\mu_i = \\mu_i^0 + RT \\ln(f_iX_i)",
"Variables": ["μi^0", "R", "T", "fi", "Xi"],
"Additional_Info": "Where μi is the chemical potential of component i, μi^0 is the standard chemical potential, R is the gas constant, T is the temperature, fi is the coefficient of chemical activity, and Xi is the molar fraction of component i."
}
]

def get_default_formulas():
    """Return the default formulas included with the package."""
    return DEFAULT_FORMULAS.copy()
