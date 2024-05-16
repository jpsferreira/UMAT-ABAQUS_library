---
title: "Title"
author: "Me"
header-includes:
   - \usepackage[bb=libus]{mathalpha}
output:
    pdf_document
---

# UMAT structure


## 1. _global_.f90
Sets the main control parameters to run the material-related routines (NELEM, NSDV, NTERM, NFACTOR).

--------------------------------------------------------------

## 2. _umat_.f90
- Main subroutine that imports _global_ parameters and organizes all secondary routines to fully describe the material behaviour.
- Header variables described in UMAT_README.txt.
- Initializes all variables as zero.

### _onem_
- Defines the 2nd and 4th order identity tensors and 4th order symmetric identity tensor: $\mathbf{I}, \mathbb{I} \mathrm{~and~} \mathbb{I}^s$, respectively. 
- $\mathbf{I}=\delta_{i j}\left(\mathbf{e}_i \otimes \mathbf{e}_j\right)$
- $\mathbb{I}=\delta_{i k} \delta_{j l}\left(\mathbf{e}_i \otimes \mathbf{e}_j \otimes \mathbf{e}_k \otimes \mathbf{e}_l\right)$
- $\mathbb{I}^s=\frac{1}{2}\left(\delta_{i k} \delta_{j l}+\delta_{i l} \delta_{j k}\right)\left(\mathbf{e}_i \otimes \mathbf{e}_j \otimes \mathbf{e}_k \otimes \mathbf{e}_l\right)$

**Material properties**: from array PROPS to scalar variable
1. k - Penalty parameter
2. c10 - Hyperelastic isotropic matrix constant 1
3. c01 - Hyperelastic isotropic matrix constant 2
4. phi - Filament volume fraction
5. ll - Filament contour length [$\mu \mathrm{m}$]
6. r0f - Filament initial end-to-end distance [$\mu \mathrm{m}$]
7. r0c - Crosslinker (CL) initial end-to-end distance [$\mu\mathrm{m}$]
8. eta - CL relative stiffness
9. mu0 - Filament stretch modulus [pN]
10. beta - Beta parameter from Holzapfel beta model
11. b0 - Bending stiffness ($T \cdot L_p \cdot k_b$) [pN * $\mu \mathrm{m}^2$]
12. lambda0 - Pre-stretch
13. nn - Isotropic filaments per unit volume [$\mu \mathrm{m}^{-3}$]
14. bb - Dispersion parameter

### _initialize_ and _sdvread_
- Initializes or reads state variables (determinant and contraction variance)

---------------------------------------------------------------

### 2.1. KINEMATICS

### _fslip_
- Computes volume-preserving/distortional part of the deformation gradient:
1. $J = \mathrm{det}(\mathbf{F})$
2. $\overline{\mathbf{F}} = J^{-1/3} \mathbf{F}$

### _matinv3d_
- This utility routine is used for computing $\mathbf{F}^{-1}$ and $\overline{\mathbf{F}}^{-1}$.

### _deformation_
- Computes right ($\mathbf{C}$) and left ($\mathbf{b}$) Cauchy-Green strain tensors, and their distortion parts.
- $\mathbf{C} = \mathbf{F}^\mathrm{T} \mathbf{F}$
- $\mathbf{b} = \mathbf{F} \mathbf{F}^\mathrm{T}$
- $\overline{\mathbf{C}} = \overline{\mathbf{F}}^\mathrm{T} \overline{\mathbf{F}}$
- $\overline{\mathbf{b}} = \overline{\mathbf{F}} \overline{\mathbf{F}}^\mathrm{T}$

### _invariants_
- Gets 1st and 2nd strain invariants ($\bar{I_i}$) [see Holzapfel 2000, 6.107-111]
- $\bar{I_1} = \mathrm{tr}\overline{\mathbf{C}} (= \mathrm{tr}\overline{\mathbf{b}})$
- $\bar{I_2} = \frac{1}{2}[(\mathrm{tr}\overline{\mathbf{C}})^2 - \mathrm{tr}(\overline{\mathbf{C}}^2)]$

### _stretch_
- Computes (deviatoric) stretch tensors $\overline{\mathbf{U}} \mathrm{~and~} \overline{\mathbf{v}}$. These tensors measure local stretching along their mutually orthogonal eigenvectors. Measure of local shape change.
- $\mathbf{U}$: Right (material) stretch tensor. Acts on the reference configuration. Its eigenvalues are the principal stretches. Shares eigenvectors ($\mathbf{\hat{N}}_a$) with $\mathbf{C}$.
- $\mathbf{v}$: Left (spatial) stretch tensor. Acts on the current configuration. Its eigenvalues are the principal stretches. Shares eigenvectors ($\mathbf{\hat{n}}_a$) with $\mathbf{b}$.
1. CALL _spectral_: ($\overline{\mathbf{C}}$) -> [$\Omega_a$ (eigenvalues = $\lambda_a^2$), $\mathbf{\hat{N}}_a$ (eigenvectors)]
2. $\bar{\lambda}_a = \sqrt{\Omega_a}$
3. $\overline{\mathbf{U}}$ in the principal referential: $\overline{\mathbf{U}}_{aa}^{\mathrm{p}} = \bar{\lambda}_a$
4. $\overline{\mathbf{U}} = \mathbf{\hat{N}} \overline{\mathbf{U}}^{\mathrm{p}} \mathbf{\hat{N}}^{\mathrm{T}}$
- Repeat, using $\overline{\mathbf{b}}$ instead of $\overline{\mathbf{C}}$, to obtain $\overline{\mathbf{v}}$.

### _rotation_
- Calculates the rotation tensor: $\mathbf{R} = \overline{\mathbf{F}} \overline{\mathbf{U}}^{-1}$

------------------------------------------------------------------------------------------

### 2.2. CONSTITUTIVE RELATIONS

### _projeul_ and _projlag_
1. Projection tensor in Lagrangian (material) description: $\mathbb{P} = \mathbb{I} - (\mathbf{C}^{-1} \otimes \mathbf{C})/3$ (Falta dividir por 3 nos pngs do UMAT-Abaqus???)
2. Projection tensor in Eulerian (spatial) description: $\mathbb{p} = \mathbb{I}^s - (\mathbf{I} \otimes \mathbf{I})/3$

### _vol_
- Volumetric contribution of the strain energy function
1. $\mathcal{G}=\frac{1}{4}\left(J^2-1-2 \ln J\right)$
2. $\Psi_{\mathrm{vol}}=k \mathcal{G}$
3. $p^{*}=\frac{\mathrm{d} \Psi_{\text {vol }}(J)}{\mathrm{d} J} = \frac{1}{2}\kappa(J-J^{-1})$
4. $\tilde{p} = p^{*} + J \frac{\mathrm{d}p^*}{\mathrm{d}J}$

#### 2.2.1. Isotropic soft ground substance (IM)
- IM contribution is only considered when filament volume fration $\varphi<1$.

##### _isomat_
- IM strain energy and its derivatives (diso) in relation to the 5 invariants.
- $\Psi_{iso}=c_{10}(\bar{I_1}-3)+c_{01}(\bar{I_2}-3)$ (Mooney-Rivlin)
- diso $= \{c_{10}, c_{01}, 0, 0, 0\}$ 

##### _pk2isomatfic_
- 2nd Piola-Kirchoff fictitious stress tensor.
1. $\bar{\gamma}_1=2\left(\frac{\partial \Psi_{\text {iso }}\left(\bar{I}_1 ; \bar{I}_2\right)}{\partial \bar{I}_1}+\bar{I}_1 \frac{\partial \Psi_{\text {iso }}\left(\bar{I}_1 ; \bar{I}_2\right)}{\partial \bar{I}_2}\right) = 2(\mathrm{diso}(1)+\bar{I}_1\mathrm{diso}(2))$
2. $\bar{\gamma}_2=-2 \frac{\partial \Psi_{\text {iso }}\left(\bar{I}_1 ; \bar{I}_2\right)}{\partial \bar{I}_2}=-2\mathrm{diso}(2)$
3. $\tilde{\mathbf{S}}=2 \frac{\partial \Psi_{\text {iso }}\left(\bar{I}_1 ; \bar{I}_2\right)}{\partial \overline{\mathbf{C}}}=\bar{\gamma}_1 \mathbf{I}+\bar{\gamma}_2 \overline{\mathbf{C}}$

##### _sigisomatfic_
- Cauchy fictitious stress tensor (push forward of $\tilde{\mathbf{S}}$).
- $ \tilde{\boldsymbol{\sigma}}=J^{-1} \mathbf{F \tilde{S} F}^{\mathrm{T}} $

##### _cmatisomatfic_
- 4th-order fictitious elasticity tensor in material description $\tilde{\mathbb{C}}$
- TO REVIEW (6.169 - Não falta multiplicar por J**-4/3?????)

##### _csisomatfic_
- 4th-order fictitious elasticity tensor in spatial description $\tilde{\mathbb{c}}$ (push-forward of $\tilde{\mathbb{C}}$)

#### 2.2.2. Filaments network (is NA implementation supposed to be here????)

##### _erfi_
- Computes the imaginary error function.

#### 2.2.3. Affine network

##### _affclnetfic_discrete_
- Affine network with compliant crosslinkers. Returns the fictitious Cauchy $\tilde{\boldsymbol{\sigma}}$ and spatial elasticity $\tilde{\mathbb{c}}$ tensors.
- Detailed in AFFCLNETFIC_README.md

##### Adding affine, non-affine and isotropic matrix contributions, for both spatial and material descriptions
- $\tilde{\boldsymbol{\sigma}}=(1-\phi)\tilde{\boldsymbol{\sigma}}_{\mathrm{IM}} +\tilde{\boldsymbol{\sigma}}_{\mathrm{NA}}+\tilde{\boldsymbol{\sigma}}_{\mathrm{AN}}$
- $\tilde{\mathbb{c}}=(1-\phi)\tilde{\mathbb{c}}_{\mathrm{IM}} +\tilde{\mathbb{c}}_{\mathrm{NA}}+\tilde{\mathbb{c}}_{\mathrm{AN}}$
- $\tilde{\mathbf{S}} = (1-\phi)\tilde{\mathbf{S}}_{\mathrm{IM}} + \tilde{\mathbf{S}}_{\mathrm{NA}}+\tilde{\mathbf{S}}_{\mathrm{AN}}$
- $\tilde{\mathbb{C}} = (1-\phi)\tilde{\mathbb{C}}_{\mathrm{IM}} + \tilde{\mathbb{C}}_{\mathrm{NA}}+\tilde{\mathbb{C}}_{\mathrm{AN}}$

##### Strain-Energy (not computed)

### 2.3. Stress measures

#### 2.3.1. Volumetric

##### _pk2vol_
- Volumetric part of the second Piola-Kirchoff stress tensor
- $\mathbf{S}_{\mathrm{vol}} = J p^* \mathbf{C}^{-1}$ (nao falta multiplicar por J????)

##### _sigvol_
- Volumetric part of the Cauchy stress tensor
- $\mathbf{\sigma}_{\mathrm{vol}} = p^* \mathbf{I}$

#### 2.3.2. Isochoric

##### _pk2iso_
- Isochoric part of the second Piola-Kirchoff stress tensor
- $\overline{\mathbf{S}} = J^{-2/3} \mathbb{P} : \tilde{\mathbf{S}} $

##### _sigiso_
- Isochoric part of the Cauchy stress tensor ($\mathbb{p}$ is the eulerian projection tensor)
- $\overline{\mathbf{\sigma}} = \mathbb{p} : \tilde{\mathbf{\sigma}}$

##### Adding volumetric and isochoric contributions
- $\mathbf{S} = \mathbf{S}_{\mathrm{vol}} + \overline{\mathbf{S}}$
- $\mathbf{\sigma} = \mathbf{\sigma}_{\mathrm{vol}} + \overline{\mathbf{\sigma}}$

### 2.4. Elasticity tensors
- Computed in both spatial and material descriptions

#### Material (not implemented)

#### Spatial

##### _setvol_
- $\mathbb{c}_{\mathrm{vol}} = \tilde{p} \mathbf{I}\otimes\mathbf{I} - 2p^* \mathbb{I}^{\mathrm{s}}$

##### _setiso_
- $\bar{\mathbb{c}} = \mathbb{p}:\tilde{\mathbb{c}}:\mathbb{p} + \frac{2}{3}tr\tilde{\mathbf{\sigma}} - \frac{2}{3}(\overline{\mathbf{\sigma}}\otimes\mathbf{I} + \mathbf{I} \otimes \overline{\mathbf{\sigma}})$

##### _setjr_
- Computes the Jaumann rate
- $\mathbb{c}^{\mathrm{jr}} = \frac{1}{2}(\mathbf{\sigma}\otimes\mathbf{I} + \mathbf{I}\otimes\mathbf{\sigma})$

##### Adding all contributions
- $\mathbb{c} = \mathbb{c}_{\mathrm{vol}} + \bar{\mathbb{c}} + \mathbb{c}^{\mathrm{jr}}$