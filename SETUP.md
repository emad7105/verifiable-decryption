# Setup
This setup is only intended for hosts with Intel or x86 architecture. It is assumed that you have a decent C compiler installed on your machine.

```bash
sudo apt-get update
sudo apt-get install libmpfr-dev
sudo apt-get install libgmp-dev
pip3 show cffi
```

This setup is performant if your CPU supports AVX512. To check this, you can run the following command.
If you find AVX512 or any mention of it in the output of the second command in the `Flags` section, your CPU supports it.

```bash
lscpu
grep -o 'avx512[^ ]*' /proc/cpuinfo
```

