// tools/vk_probe.cpp
// 只读探测本机 Vulkan 设备能力，用于确定 compute 后端的实现方案。
// 构建（VS 开发者环境）：
//   cl /nologo /EHsc /std:c++17 /I"%VULKAN_SDK%\Include" tools\vk_probe.cpp ^
//      /Fe:tools\vk_probe.exe /link /LIBPATH:"%VULKAN_SDK%\Lib" vulkan-1.lib
#include <vulkan/vulkan.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#define PR(cond) std::printf("  %-46s : %s\n", #cond, (cond) ? "YES" : "no")

static std::vector<VkExtensionProperties> deviceExtensions(VkPhysicalDevice d) {
  uint32_t n = 0; vkEnumerateDeviceExtensionProperties(d, nullptr, &n, nullptr);
  std::vector<VkExtensionProperties> v(n);
  vkEnumerateDeviceExtensionProperties(d, nullptr, &n, v.data());
  return v;
}
static bool has(const std::vector<VkExtensionProperties>& v, const char* name) {
  for (auto& e : v) if (!std::strcmp(e.extensionName, name)) return true;
  return false;
}

int main() {
  uint32_t api = VK_API_VERSION_1_3;
  VkApplicationInfo ai{VK_STRUCTURE_TYPE_APPLICATION_INFO};
  ai.pApplicationName = "vk_probe"; ai.apiVersion = api;
  VkInstanceCreateInfo ci{VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
  ci.pApplicationInfo = &ai;
  VkInstance inst = VK_NULL_HANDLE;
  if (vkCreateInstance(&ci, nullptr, &inst) != VK_SUCCESS) { std::printf("instance failed\n"); return 1; }

  std::printf("=== instance ===\n");
  uint32_t iv = 0; vkEnumerateInstanceVersion(&iv);
  std::printf("  loader version                     : %u.%u.%u\n",
              VK_VERSION_MAJOR(iv), VK_VERSION_MINOR(iv), VK_VERSION_PATCH(iv));

  uint32_t nd = 0; vkEnumeratePhysicalDevices(inst, &nd, nullptr);
  std::vector<VkPhysicalDevice> devs(nd);
  vkEnumeratePhysicalDevices(inst, &nd, devs.data());
  std::printf("  physical devices                   : %u\n\n", nd);

  for (uint32_t i = 0; i < nd; ++i) {
    VkPhysicalDevice d = devs[i];
    VkPhysicalDeviceProperties p{};
    vkGetPhysicalDeviceProperties(d, &p);
    VkPhysicalDeviceShaderFloat16Int8Features f16{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SHADER_FLOAT16_INT8_FEATURES};
    VkPhysicalDevice16BitStorageFeatures f16s{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_16BIT_STORAGE_FEATURES, &f16};
    VkPhysicalDevice8BitStorageFeatures f8s{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_8BIT_STORAGE_FEATURES, &f16s};
    VkPhysicalDeviceVulkan12Features f12{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_2_FEATURES, &f8s};
    VkPhysicalDeviceVulkan11Features f11{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_VULKAN_1_1_FEATURES, &f12};
    VkPhysicalDeviceFeatures2 f2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_FEATURES_2, &f11};
    vkGetPhysicalDeviceFeatures2(d, &f2);

    VkPhysicalDeviceSubgroupProperties sg{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_SUBGROUP_PROPERTIES};
    VkPhysicalDeviceProperties2 p2{VK_STRUCTURE_TYPE_PHYSICAL_DEVICE_PROPERTIES_2, &sg};
    vkGetPhysicalDeviceProperties2(d, &p2);

    std::printf("=== device %u: %s ===\n", i, p.deviceName);
    std::printf("  api %u.%u.%u | driver %u | type %d | subgroupSize %u\n",
                VK_VERSION_MAJOR(p.apiVersion), VK_VERSION_MINOR(p.apiVersion), VK_VERSION_PATCH(p.apiVersion),
                p.driverVersion, (int)p.deviceType, sg.subgroupSize);

    std::printf("--- features (对本项目关键) ---\n");
    PR(f2.features.shaderFloat64);
    PR(f2.features.shaderInt64);
    PR(f16.shaderFloat16);
    PR(f16.shaderInt8);
    PR(f16s.storageBuffer16BitAccess);
    PR(f16s.uniformAndStorageBuffer16BitAccess);
    PR(f8s.storageBuffer8BitAccess);
    PR(f12.shaderSharedInt64Atomics);
    PR(f12.shaderBufferInt64Atomics);
    PR(f12.bufferDeviceAddress);
    PR(f12.timelineSemaphore);
    PR(f12.hostQueryReset);
    PR(f11.storageBuffer16BitAccess);
    PR(f11.variablePointersStorageBuffer);

    std::printf("--- subgroup ---\n");
    PR(sg.supportedStages & VK_SHADER_STAGE_COMPUTE_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_BASIC_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_VOTE_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_ARITHMETIC_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_BALLOT_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_SHUFFLE_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_SHUFFLE_RELATIVE_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_CLUSTERED_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_QUAD_BIT);
    PR(sg.supportedOperations & VK_SUBGROUP_FEATURE_PARTITIONED_BIT_NV);
    std::printf("  quadOperationsInAllStages          : %s\n", sg.quadOperationsInAllStages ? "YES" : "no");
    std::printf("  subgroupSupportedOperations mask   : 0x%x\n", (unsigned)sg.supportedOperations);

    std::printf("--- limits (compute) ---\n");
    std::printf("  maxComputeSharedMemorySize         : %u B (%.1f KB)\n",
                p.limits.maxComputeSharedMemorySize, p.limits.maxComputeSharedMemorySize / 1024.0);
    std::printf("  maxComputeWorkGroupInvocations     : %u\n", p.limits.maxComputeWorkGroupInvocations);
    std::printf("  maxComputeWorkGroupSize            : %u %u %u\n",
                p.limits.maxComputeWorkGroupSize[0], p.limits.maxComputeWorkGroupSize[1], p.limits.maxComputeWorkGroupSize[2]);
    std::printf("  maxComputeWorkGroupCount           : %u %u %u\n",
                p.limits.maxComputeWorkGroupCount[0], p.limits.maxComputeWorkGroupCount[1], p.limits.maxComputeWorkGroupCount[2]);
    std::printf("  maxStorageBufferRange              : %llu MB\n",
                (unsigned long long)(p.limits.maxStorageBufferRange >> 20));
    std::printf("  maxDescriptorSetStorageBuffers     : %u\n", p.limits.maxDescriptorSetStorageBuffers);
    std::printf("  maxDescriptorSetStorageBuffersDyn  : %u\n", p.limits.maxDescriptorSetStorageBuffersDynamic);
    std::printf("  maxPerStageDescriptorStorageBufs   : %u\n", p.limits.maxPerStageDescriptorStorageBuffers);
    std::printf("  maxBoundDescriptorSets             : %u\n", p.limits.maxBoundDescriptorSets);
    std::printf("  maxPushConstantsSize               : %u B\n", p.limits.maxPushConstantsSize);
    std::printf("  maxMemoryAllocationCount           : %u\n", p.limits.maxMemoryAllocationCount);

    std::printf("--- memory heaps ---\n");
    VkPhysicalDeviceMemoryProperties mp{};
    vkGetPhysicalDeviceMemoryProperties(d, &mp);
    for (uint32_t h = 0; h < mp.memoryHeapCount; ++h)
      std::printf("  heap %u: %.2f GB %s\n", h, mp.memoryHeaps[h].size / 1073741824.0,
                  (mp.memoryHeaps[h].flags & VK_MEMORY_HEAP_DEVICE_LOCAL_BIT) ? "(device local)" : "(host visible)");

    auto ext = deviceExtensions(d);
    std::printf("--- extensions (%zu total, 关键项) ---\n", ext.size());
    const char* k[] = {
      "VK_EXT_shader_atomic_float", "VK_EXT_shader_atomic_float2",
      "VK_KHR_shader_atomic_int64", "VK_EXT_shader_subgroup_ballot",
      "VK_KHR_shader_float16_int8", "VK_KHR_16bit_storage",
      "VK_KHR_cooperative_matrix", "VK_NV_cooperative_matrix",
      "VK_KHR_synchronization2", "VK_KHR_timeline_semaphore",
      "VK_EXT_host_query_reset", "VK_KHR_buffer_device_address",
      "VK_KHR_shader_clock", "VK_KHR_shader_integer_dot_product",
      "VK_EXT_memory_budget", "VK_EXT_memory_priority",
      "VK_KHR_external_memory_win32", "VK_EXT_external_memory_host",
      "VK_KHR_external_semaphore_win32", "VK_EXT_descriptor_indexing",
      "VK_KHR_maintenance4", "VK_KHR_maintenance5", "VK_KHR_maintenance6",
      "VK_KHR_dynamic_rendering",
    };
    for (auto* e : k) std::printf("  %-46s : %s\n", e, has(ext, e) ? "YES" : "no");
    std::printf("\n");
  }
  vkDestroyInstance(inst, nullptr);
  return 0;
}
