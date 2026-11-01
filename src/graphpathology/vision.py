"""Explicit frozen encoder construction. Downloads are opt-in."""
def resnet18(*, download=False):
    from torchvision.models import resnet18 as make, ResNet18_Weights
    from torch import nn
    model=make(weights=ResNet18_Weights.DEFAULT if download else None)
    model.fc=nn.Identity()
    model.requires_grad_(False)
    return model.eval()

def patch_transform(image):
    from PIL import Image
    from torchvision.models import ResNet18_Weights
    image=image.convert('RGB').resize((224,224),Image.Resampling.LANCZOS)
    return ResNet18_Weights.DEFAULT.transforms()(image)

def ensemble(heads, features, feature_scaler, target_scaler):
    import torch
    import numpy as np
    if len(heads)!=3: raise ValueError('Exactly three original heads required')
    x=torch.from_numpy(feature_scaler.transform(features))
    with torch.inference_mode():
        return np.mean([target_scaler.inverse(h.eval()(x).numpy()) for h in heads],axis=0)
