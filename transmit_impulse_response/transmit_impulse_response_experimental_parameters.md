# Transmit Impulse Response: Experimental Parameters, Data Extraction and Calibration

Hydrophone linescan measurements to measure the transmit impulse response of individual open-UST transducer elements.

## Data Access
Data folder: `Datasets\open-UST-imaging-and-evaluation\transmit_impulse_response\`

Raw voltage data subfolders:
|Sub-foldername|Description|
|----------|-------------|
|`\Impulse Response 80ns 90V 220mm Linescans Probe Q`|module Q (without acoustic matching layers)|
|`\Impulse Response 80ns 90V 220mm Linescans with 22uH matching`|module A driven thorugh electrical impedance matching inductors|
|`\Impulse Response 80ns 90V 220mm Linescans`|modules A, E, F and G (with acoustic matching layers)|

Calibrated pressure data filenames: 
|Filename|Description|
|----------|-------------|
|`transmit_impulse_response_probe_q.mat`|module Q (without acoustic matching layers)|
|`transmit_impulse_response_probes_A22uH.mat`|module A driven thorugh electrical impedance matching inductors|
|`transmit_impulse_response_probes_AEFG.mat`|modules A, E, F and G (with acoustic matching layers)|

---
## Setup
### Measurement Setup
- Tank of deionised water
- X-Y-Z stepper motor controlled positioning system (Precision Acoustics)

### Receive signal chain
1. 0.2mm needle hydrophone (Precision Acoustics 3492)
1. Submersible preamplifier
1. DC Coupler (Precision Acoustics DCPS340)
1. 50 ohm shunt
1. Digital oscilloscope Channel 2 (Tektronix, DPO5034B)

### Scope settings
- dt = 2 ns/pt
- Fs = 500 MS/s
- 500 averages
- Falling-edge trigger into channel 1

### Transmit settings
- Verasonics Vantage 256
- `TPC(1).hv = 90` (90V drive voltage)
- `TW(1).PulseCode = [0, 20, 0, 0, 1]` (pulse width of 80 ns, or 20 clock cycles at 250 MHz)
- `SOURCE_PULSE_REPETITION = 0.001` (1 ms source repetition period)

### Acoustic settings
- Time-of-flight at beam-max = 148 us
- z-position at beam max = 220 mm
- Temperature = 22.1 degC (measured with NI J-type thermocouple)

---
## Measurement Procedure
Interelement variation in the transmit-impulse-response is to be measured for a range of transducer elements. To allow direct comparison, the impulse response waveform must be measured at the same point in the field for each transducer element. Due to skew in the elevation and lateral planes, the hydrophone must be well aligned with the beam axis for each element. It is assumed that:
- The face of the transducer module is coplanar with the X-Y axes of the positioning system.
- The beam skew is small enough not to significantly change the acoustic path length (which would lead to an increase in geometric spreading and an underestimate of pressure).
- Therefore, the z-position does not need to be corrected for each element.

For each element a linescan is acquired so that the optimal x-position can be found later (the UMS software auto-beam-alignment can fail sometimes).

Setup:
- Mount the transducer module in the tank of deionised water.
- Mount the needle hydrophone to the arm of the positioning system.
- Align the hydrophone with the centre of the transducer module, visually.
- Transmit on channel 1 of the module.
- Use the time of flight on the scope to estimate the z-position of the hydrophone.
- Move the hydrophone to the required z-position.

Align to the first element in each module:
- Adjust the Verasonics transmit settings so the correct channel is turned on.
- Roughly align the hydrophone with the beam axis:
    - Peak value of voltage squared integral
    - N = 60 pts
    - dx = 0.25 mm
    - X/Y/X, 1.5 iterations
- Set the acquisition time window of the scope

Then, for each element:
- Adjust the Verasonics transmit settings so the correct channel is turned on.
- Move in the x-direction by the element pitch (2.54 mm).
- Relocate the maximum along the y-axis:
    - Peak value of voltage squared integral
    - N = 10 pts
    - dx = 0.25 mm
    - Y/Y, 1 iteration
- Acquire data: collect a linescan in the x-direction (lateral plane):
    - N = 160 pts
    - dx = 0.25 mm

---
## Data Extraction and Calibration

The measured voltage is converted to pressure by deconvolving the frequency response of the hydrophone (using sensitivity magnitude data). A boxcar filter is first applied to the voltage data with width equal to the calibration bandwidth.

Extraction and calibration script: `\extract_calibrate_transmit_impulse_response.m`