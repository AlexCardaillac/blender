/* SPDX-FileCopyrightText: 2005 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "node_shader_util.hh"

#include "BLI_string.h"

#include "UI_interface.hh"
#include "UI_resources.hh"

#include "BKE_node_runtime.hh"

namespace blender::nodes::node_shader_volume_absorption_cc {

static void node_declare(NodeDeclarationBuilder &b)
{
  b.add_input<decl::Color>("Color").default_value({0.8f, 0.8f, 0.8f, 1.0f});
#define SOCK_COLOR_ID 0
  b.add_input<decl::Float>("Density").default_value(1.0f).min(0.0f).max(1000.0f);
#define SOCK_DENSITY_ID 1
  b.add_input<decl::Vector>("Densities")
      .default_value({1.0f, 1.0f, 1.0f})
      .min(0.0f)
      .max(1000.0f)
      .subtype(PROP_FACTOR);
  b.add_input<decl::Float>("Weight").unavailable();
  b.add_output<decl::Shader>("Volume").translation_context(BLT_I18NCONTEXT_ID_ID);
}

static void node_shader_buts_absorption(uiLayout *layout, bContext * /*C*/, PointerRNA *ptr)
{
  uiItemR(layout, ptr, "density_mode", UI_ITEM_R_SPLIT_EMPTY_NAME, "", ICON_NONE);
}

static void node_shader_init_absorption(bNodeTree * /*ntree*/, bNode *node)
{
  node->custom1 = SHD_DENSITY_GLOBAL;
}

static void node_shader_update_absorption(bNodeTree *ntree, bNode *node)
{
  const int density_mode = node->custom1;

  LISTBASE_FOREACH (bNodeSocket *, sock, &node->inputs) {
    if (STR_ELEM(sock->name, "Density")) {
      bke::nodeSetSocketAvailability(ntree, sock, density_mode == SHD_DENSITY_GLOBAL);
    }
    else if (STR_ELEM(sock->name, "Color")) {
      bke::nodeSetSocketAvailability(ntree, sock, density_mode == SHD_DENSITY_GLOBAL);
    }
    else if (STR_ELEM(sock->name, "Densities")) {
      bke::nodeSetSocketAvailability(ntree, sock, density_mode == SHD_DENSITY_CHANNEL);
    }
  }
}

static int node_shader_gpu_volume_absorption(GPUMaterial *mat,
                                             bNode *node,
                                             bNodeExecData * /*execdata*/,
                                             GPUNodeStack *in,
                                             GPUNodeStack *out)
{
  if (node_socket_not_zero(in[SOCK_DENSITY_ID]) && node_socket_not_white(in[SOCK_COLOR_ID])) {
    GPU_material_flag_set(mat, GPU_MATFLAG_VOLUME_ABSORPTION);
  }
  return GPU_stack_link(mat, node, "node_volume_absorption", in, out);
}

#undef SOCK_COLOR_ID
#undef SOCK_DENSITY_ID

}  // namespace blender::nodes::node_shader_volume_absorption_cc

/* node type definition */
void register_node_type_sh_volume_absorption()
{
  namespace file_ns = blender::nodes::node_shader_volume_absorption_cc;

  static bNodeType ntype;

  sh_node_type_base(&ntype, SH_NODE_VOLUME_ABSORPTION, "Volume Absorption", NODE_CLASS_SHADER);
  ntype.declare = file_ns::node_declare;
  ntype.add_ui_poll = object_shader_nodes_poll;
  ntype.draw_buttons = file_ns::node_shader_buts_absorption;
  ntype.initfunc = file_ns::node_shader_init_absorption;
  ntype.updatefunc = file_ns::node_shader_update_absorption;
  ntype.gpu_fn = file_ns::node_shader_gpu_volume_absorption;
  blender::bke::node_type_size_preset(&ntype, blender::bke::eNodeSizePreset::MIDDLE);

  nodeRegisterType(&ntype);
}
