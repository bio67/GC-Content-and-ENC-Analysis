import pandas as pd

# Read input file
input_file = "codonw_results.csv"  # Replace with your file name
data = pd.read_csv(input_file)

# Define a function to calculate ENC_exp
def calculate_ENC_exp(gc3s):
    # Avoid division by zero when GC3s is 0 or 1
    if gc3s <= 0:
        gc3s = 0.0001
    elif gc3s >= 1:
        gc3s = 0.9999
    return 2 + (29 / (gc3s**2 + (1 - gc3s)**2))

# Calculate ENC_exp and (ENC_exp - ENC_obs) / ENC_exp
data['ENC_exp'] = data['GC3s'].apply(calculate_ENC_exp)
data['(ENC_exp - ENC_obs) / ENC_exp'] = (data['ENC_exp'] - data['ENC_obs']) / data['ENC_exp']

# Save results to a new file
output_file = "ENCratio_result.csv"
data.to_csv(output_file, index=False)

print(f"Results saved to {output_file}")