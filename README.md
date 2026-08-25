# EKG Maker 🫀

**EKG Maker** is a web-based interactive Electrocardiogram (ECG/EKG) simulator built for medical education and clinical simulation. Designed by a physician, it generates real-time, customizable 12-lead EKG tracings to help students and professionals master cardiovascular diagnostics.

[![Live Demo / Website](https://img.shields.io/badge/Demo-Coming_Soon-blue.svg)](#)
[![Tech Stack: Vanilla JS / HTML / CSS](https://img.shields.io/badge/Tech-Vanilla_JS-yellow.svg)](#)

---

## 🚀 Overview

**Why I built this:**  
As a medical doctor, I faced the challenge of providing students and colleagues with accessible, realistic, and dynamic cardiology clinical cases. Traditional static images are limited. I taught myself to code and built **EKG Maker** to solve this problem programmatically.

**What this project demonstrates:**
- **Problem-Solving & Domain Expertise:** Bridging the gap between complex medical requirements (cardiology electrophysiology) and software development.
- **Algorithms & Math:** Implemented procedural generation for EKG waveforms using mathematical modeling, Bezier-like curves, and interpolation for realistic heartbeat rendering.
- **Frontend Architecture:** Built from scratch using **Vanilla JavaScript, HTML5 Canvas, and CSS3** (No heavy frameworks). It demonstrates a strong grasp of core web technologies, state management, and real-time rendering.
- **Performance & Optimization:** Optimized Canvas API drawing for a smooth real-time "monitor" view that runs at 60 FPS.
- **UI/UX Design:** Designed a clean, intuitive, and responsive dashboard that presents dense medical parameters in an accessible way.

---

## ✨ Key Features

- **12-Lead Simulation:** Generates accurate 12-lead EKG views (I, II, III, aVR, aVL, aVF, V1-V6).
- **Real-Time Monitor:** Features a live, animated EKG monitor with adjustable speed (25/50 mm/s) and gain.
- **Pathology Generation:** Customize rhythms, heart axes, and waveforms to simulate Left/Right Bundle Branch Blocks (LBBB/RBBB), Ischemia (ST elevation/depression), and artifacts.
- **Export & Print:** Export generated cases as JSON files or print them in standard medical EKG grid formats (Letter/A4).
- **Responsive Layout:** Adaptive CSS Grid and Flexbox layout for various screen sizes.

---

## 🛠 Tech Stack

- **Frontend:** HTML5, CSS3 (Variables, Grid, Flexbox)
- **Logic & State:** Vanilla JavaScript (ES6+)
- **Rendering:** HTML5 Canvas API

---

## 📂 Project Structure

- `index.html`: Main layout and UI structure.
- `style.css`: UI styling, dark mode theme, and printable formats.
- `app.js`: Core engine containing wave generation math, state management, DOM manipulation, and Canvas rendering.

---

## 💻 How to Run Locally

Since this is a client-side only application without build steps or external dependencies, running it is incredibly simple:

1. Clone the repository:
   ```bash
   git clone https://github.com/arzaluz-chris/EKG-Maker.git
   ```
2. Navigate to the directory:
   ```bash
   cd EKG-Maker
   ```
3. Open `index.html` in your favorite browser. No local server required!

---

## 🤝 Let's Connect

I am a physician transitioning into tech, passionate about HealthTech, software engineering, and building tools that make a difference. 

- **LinkedIn:** [My LinkedIn](https://www.linkedin.com/in/christian-arzaluz/)
- **GitHub:** [arzaluz-chris](https://github.com/arzaluz-chris)

*Feedback, suggestions, and networking are always welcome!*
