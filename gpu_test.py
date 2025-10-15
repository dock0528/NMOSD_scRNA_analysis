import torch

def test_gpu():
    
    print("PyTorch version:", torch.__version__)
    is_cuda_available = torch.cuda.is_available()
    print("CUDA Available:", is_cuda_available)
    
    if is_cuda_available:

        print("GPU Name:", torch.cuda.get_device_name(0))
        device = torch.device("cuda") 
        print("Running on Device:", device)
        tensor = torch.rand((1000, 1000)).to(device)
        print("Tensor created on GPU")
        result = tensor @ tensor
        print("Matrix multiplication successful on GPU")

    else:
        print("No GPU detected. Please check your CUDA setup.")

if __name__ == "__main__":
    test_gpu()
