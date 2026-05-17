# Git Commands to Initialize and Push Repository

Run the following commands in a terminal from the repository root:

```bash
cd github_release/CCA_CAF_TAM_Communication

# Initialize git
git init

# Stage all files
git add .

# First commit
git commit -m "Initial release of analysis scripts for CCA CAF-TAM communication study"

# Set main branch
git branch -M main

# Add remote (REPLACE YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/CCA_CAF_TAM_Communication.git

# Push
git push -u origin main
```

## After Pushing

1. Verify the repository is public on GitHub.
2. Generate a Zenodo DOI for the repository snapshot.
3. Update the manuscript Data and Code Availability section:

```
All analysis scripts are available at https://github.com/YOUR_USERNAME/CCA_CAF_TAM_Communication
(Zenodo DOI: 10.5281/zenodo.XXXXXXX).
```

4. Replace `YOUR_USERNAME` and the Zenodo DOI with actual values before manuscript submission.
