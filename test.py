import torch
import time

def test_gpu_with_torch():
    print("===== PyTorch GPU 測試 =====")
    if not torch.cuda.is_available():
        print("❌ 沒有偵測到 GPU")
        return
    
    device = torch.device("cuda")
    print(f"✅ 偵測到 GPU: {torch.cuda.get_device_name(device)}")
    
    # 建立隨機矩陣
    a = torch.randn(10000, 10000, device=device)
    b = torch.randn(10000, 10000, device=device)

    torch.cuda.synchronize()
    start = time.time()

    # 矩陣相乘
    c = torch.matmul(a, b)

    torch.cuda.synchronize()
    end = time.time()

    print(f"運算完成，耗時: {end - start:.4f} 秒")

if __name__ == "__main__":
    test_gpu_with_torch()
