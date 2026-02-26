Transforms from eer to mrc with super-resolution, motion corrected.
Designed by Shahar Seifer, Elbaum lab, Weizmann Institute of Science, 2026
Based on https://github.com/fei-company/EerReaderLib/tree/master

Usage: Either run the exe file from a folder that contains "input.eer" or use the full structure:
eer2mrc [0] [input_file.eer] [output_file.mrc] [upscaleFactor] [timestep_S] [pixelSize_A]
For debugging:
eer2mrc 1