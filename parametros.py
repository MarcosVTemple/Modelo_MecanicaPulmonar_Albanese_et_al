cts_mp = {
    "f": 0.20,  # [rad/s]
    "N": 50000,
    "dt": 0.0001,

    # COMPLIANCE l/cmH2O
    "Cl": 0.00127,
    "Ctr": 0.00238,
    "Ccw": 0.2445,
    "Cb": 0.0131,
    "CA": 0.2,
    # RESISTANCE cmH2O·s/l
    "Rtb": 0.3063,
    "Rlt": 0.3369,
    "Rml": 1.021,
    "RbA": 0.0817,
    # UNSTRESSED VOLUME L
    "Vul": 0.0344,
    "Vut": 0.00663,
    "Vub": 0.0187,
    "VuA": 1.263,
    # ADDITIONAL PARAMETERS
    "IEratio": 0.6,
    "RR": 12,   # breaths/min
    "Pmus_min": -5,  # cmH2O
    "Ppl_EE": -5,  # cmH2O

    "Ti": 1.875,  # [s]
    "Te": 3.125,  # [s]
    "T": 5,   # [s]
    # "tau": Te/5,
    "Pao": 0,  #At the beginning of inspiration, alveolar pressure equals Patm, i.e., zero pressure,  # cmH2O
    "Pvent": 0,
}