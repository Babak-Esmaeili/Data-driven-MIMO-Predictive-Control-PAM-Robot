
# Data-Driven MIMO Predictive Model-Free Adaptive Integral Terminal Sliding Mode Control for PAM-Driven Robotic Manipulators

This repository contains MATLAB code and simulation setups for our paper:

📄 **Title**: Data-driven MIMO discrete-time predictive model-free adaptive integral terminal sliding mode controller design for robotic manipulators driven by pneumatic artificial muscles  
📰 **Conference**: 2019 6th International Conference on Control, Instrumentation and Automation (ICCIA)    
🔗 [DOI: 10.1109/ICCIA49288.2019.9030938](https://doi.org/10.1109/ICCIA49288.2019.9030938)

---

## 🧠 Abstract

In this paper, a data-driven multi-input multi-output discrete-time predictive model-free adaptive integral terminal sliding mode controller is proposed for the robotic manipulators driven by pneumatic artificial muscles with unknown models and external disturbances. First, using the concept of compact-form dynamic linearization technique, dynamics of the robotic manipulator are represented in discrete-time linear data-model form. Second, a time-delay based estimation algorithm along with an adaptive law is investigated to estimate external disturbances and pseudo-Jacobian matrix, respectively. Eventually, the proposed controller is derived based on a data-driven discrete-time integral terminal sliding function and by utilizing the concept of model-predictive controllers. The simulation results obviously clarify the efficacy of the proposed method.

---

## 🎯 Overview

This work proposes a robust data-driven control strategy for rigid robotic manipulators actuated by pneumatic artificial muscles (PAMs) with unknown dynamics and disturbances. Key contributions include:
- **MIMO discrete-time predictive controller** using model-free adaptive integral terminal sliding mode concepts.
- **Dynamic linearization** to derive a linear data-model form via pseudo-Jacobian estimation with adaptive law.
- **Disturbance estimation** using a one-step delayed algorithm.
- **Finite-time convergence** ensured by nonlinear integral terminal sliding surfaces.
- **Predictive component** via receding-horizon optimization for improved tracking accuracy and robustness.

---

## 🛠 Requirements

- MATLAB R2019b or newer  
- No additional toolboxes required  

---

## ▶️ Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/Babak-Esmaeili/Data-Driven-MIMO-Predictive-Control-PAM-Robot.git
   cd Data-Driven-MIMO-Predictive-Control-PAM-Robot/Codes
   ```
2. Open MATLAB and navigate to the `Codes/` folder.  
3. Run `main.m` to reproduce:
   - Joint tracking performance  
   - Control input smoothness and variation  
   - Sliding surface evolution and predictive improvements

---

## 📜 License and Contact

This project is licensed under the MIT License – see [LICENSE](LICENSE).  
For questions or collaboration, contact:

- **Babak Esmaeili** – esmaeil1@msu.edu

---

## 📚 Citation

If you use this repository, please cite:

```bibtex
@inproceedings{esmaeili2019mimo,
  title={Data-driven MIMO discrete-time predictive model-free adaptive integral terminal sliding mode controller design for robotic manipulators driven by pneumatic artificial muscles},
  author={Esmaeili, Babak and Baradarannia, Mahdi and Salim, Mina and Farzamnia, Ali},
  booktitle={7th RSI International Conference on Robotics and Mechatronics (ICRoM)},
  year={2019},
  organization={IEEE},
  doi={10.1109/ICRoM48714.2019.9071819}
}
```

---
