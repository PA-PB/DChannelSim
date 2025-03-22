# 📡 Digital Communication System Simulation – ELC1120

This project simulates a digital communication system, developed for the **Telecomunicações II ELC1120** course. The goal is to evaluate how digital data can be transmitted over a noisy channel and recovered with minimal errors.

---

## 🧩 The Problem

How can we **transmit a digital message reliably** through a **noisy communication channel**, ensuring the integrity of the original message?

---

## 🛠️ The Solution

The solution involves simulating a full digital communication pipeline with encoding, modulation, and error correction techniques:

1. **Information Source**  
   - Text: *Alice in Wonderland*

2. **Source Encoding**  
   - Huffman coding for efficient compression

3. **Channel Encoding**  
   - Hamming code (13,8) to enable error detection and correction

4. **Modulation**  
   - 4-PSK (Phase Shift Keying) maps bit pairs to complex symbols

5. **Noisy Channel**  
   - AWGN (Additive White Gaussian Noise) simulates real-world interference

6. **Demodulation & Decoding**  
   - 4-PSK demodulation  
   - Hamming decoding for error correction  
   - Huffman decoding to reconstruct the original message

7. **Evaluation**  
   - BER (Bit Error Rate) is measured with and without correction  
   - Visual analysis through constellation plots

---

## 📈 Results

- Error correction significantly improves performance
- At **Es/No ≥ 11 dB**, the BER reaches **zero**
- Even under noise, the system can **recover the original text**

---

## 📎 Full Report

[📄 Download the full PDF report](./RelatorioFinal_Telecom_PedroAPB.pdf)
