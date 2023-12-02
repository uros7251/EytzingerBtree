#pragma once
// -------------------------------------------------------------------------------------
// this file is taken from https://github.com/leanstore/leanstore/blob/master/backend/leanstore/storage/buffer-manager/Swip.hpp
// and modified
// -------------------------------------------------------------------------------------
#include <atomic>
#include <cassert>
// -------------------------------------------------------------------------------------
namespace guidedresearch
{
typedef uint64_t u64;
typedef uint8_t u8;
// -------------------------------------------------------------------------------------
class BufferFrame;  // Forward declaration
// -------------------------------------------------------------------------------------
class Swip
{
   // -------------------------------------------------------------------------------------
   // 1xxxxxxxxxxxx unswizzled, 0xxxxxxxxxxx swizzled
   static const u64 unswizzled_bit = u64(1) << 63;
   static const u64 mask = ~(u64(1) << 63);
   static_assert(unswizzled_bit == 0x8000000000000000, "");
   static_assert(mask == 0x7FFFFFFFFFFFFFFF, "");

   public:
   union {
      u64 pid;
      BufferFrame* bf;
   };
   // -------------------------------------------------------------------------------------
   Swip() = default;
   Swip(BufferFrame* bf) : bf(bf) {}
   Swip(u64 pid) : pid(pid | unswizzled_bit) {}
   template <typename T2>
   Swip(Swip& other) : pid(other.pid)
   {
   }
   // -------------------------------------------------------------------------------------
   bool operator==(const Swip& other) const { return (raw() == other.raw()); }
   // -------------------------------------------------------------------------------------
   bool isSWIZZLED() const { return !isUNSWIZZLED(); }
   bool isUNSWIZZLED() const { return pid & unswizzled_bit; }
   // -------------------------------------------------------------------------------------
   u64 asPageID() const { return pid & mask; }
   BufferFrame& asBufferFrame() const { return *bf; }
   BufferFrame& asBufferFrameMasked() const { return *reinterpret_cast<BufferFrame*>(pid & mask); }
   u64 raw() const { return pid; }
   // -------------------------------------------------------------------------------------
   void warm(BufferFrame* bf) { this->bf = bf; }
   // -------------------------------------------------------------------------------------
   void cool(u64 pid) { this->pid = pid | unswizzled_bit; }
   // -------------------------------------------------------------------------------------
   static Swip fromPID(u64 pid) { return Swip(pid); }
};
// -------------------------------------------------------------------------------------
}  // namespace guidedresearch